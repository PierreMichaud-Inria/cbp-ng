#include <cstdlib>
#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;

// this is a basic TAGE, not necessarily well optimized


#define USE_META
#define RESET_UBITS


template<u64 NUMG, u64 LOGG, u64 LOGB, u64 TAGW, u64 GHIST>
struct tage : predictor {
  static constexpr u64 CTR = 3;
  static constexpr u64 MINHIST = 3;
  static constexpr u64 WEAK1 = 1<<(CTR-1);
  static constexpr u64 WEAK0 = WEAK1-1;
  static constexpr u64 METABITS = 4;
  static constexpr u64 UCTRBITS = 8;
  static constexpr u64 PATHBITS = 6;
#ifdef USE_META
  static constexpr u64 METAPIPE = 2;
#endif
  static_assert(PATHBITS > loglineinst);
  static_assert(LOGB > loglineinst);
  static constexpr u64 BINDEXBITS = LOGB-loglineinst;
  static_assert(TAGW > loglineinst); // the unhashed line offset is part of the tag
  static constexpr u64 HTAGBITS = TAGW-loglineinst; // hashed tag bits

  geometric_folds<NUMG,MINHIST,GHIST,LOGG,HTAGBITS> gfolds;

  reg<64-loglineinst> lineaddr;
  reg<BINDEXBITS> bindex; // bimodal table index 
  arr<reg<LOGG>,NUMG> gindex; // global tables indexes
  arr<reg<HTAGBITS>,NUMG> htag; // computed hashed tags

  arr<reg<2>,lineinst> readb; // read bimodal counters for each offset
  arr<reg<TAGW>,NUMG> readt; // read tags
  arr<reg<CTR>,NUMG> readc; // read global counters
  arr<reg<1>,NUMG> readu; // read u bits
  reg<NUMG> notumask; // read u bits, inverted

  arr<reg<NUMG+1>,lineinst> match; // all matches for each offset
  arr<reg<NUMG+1>,lineinst> match1; // longest match for each offset
  arr<reg<NUMG+1>,lineinst> match2; // second longest match for each offset

  arr<reg<1>,lineinst> pred1; // primary prediction for each offset
  arr<reg<1>,lineinst> pred2; // alternate prediction for each offset
  arr<reg<1>,lineinst> prediction; // final prediction for each offset

#ifdef USE_META
  arr<reg<METABITS,i64>,METAPIPE> meta; // select between pred1 and pred2
  arr<reg<1>,lineinst> newly_alloc;
#endif

#ifdef RESET_UBITS
  reg<UCTRBITS> uctr; // u bits counter (reset u bits when counter saturates)
#endif

  ram<val<TAGW>,(1<<LOGG)> gtag[NUMG] {"TAGS"}; // global tables tags
  ram<val<CTR>,(1<<LOGG)> gctr[NUMG] {"3-BIT CTRS"}; // global tables counters
  ram<val<2>,(1<<BINDEXBITS)> bim[lineinst] {"2-BIT CTRS"}; // bimodal table

  zone UPDATE_ONLY;
  ram<val<1>,(1<<LOGG)> ubit[NUMG] {"U BITS"}; // "useful" bits

  // simulation artifacts (hardware cost may not be real)
  u64 update_rank = 0;
  arr<reg<loglineinst>,lineinst> update_offset; // line offset of each branch in the block
  arr<reg<1>,lineinst> update_dir; // direction of each branch in the block
  arr<reg<PATHBITS>,lineinst> update_nextinst; // next inst address of each branch in the block
  arr<reg<lineinst>,lineinst> update_mask; // for each offset, which branch of the block has that offset
  arr<reg<1>,lineinst> update_last_branch; // the 1 bit tells the offset of the last branch in the block
  reg<loglineinst> last_offset; // line offset of last branch in the block
  reg<1> last_dir; // direction of last branch in the block
  reg<PATHBITS> target_inst; // bits to inject in path history if block exits with jump

  tage()
  {
    if (LOGG == HTAGBITS) {
      std::cout << "WARNING: the tag function and index function are not different enough\n";
    }
    std::cout << "TAGE global entries: " << NUMG << " x " << (1<<LOGG) << std::endl;
    std::cout << "TAGE bimodal entries: " << (1<<LOGB) << std::endl;
    std::cout << "TAGE global entry tag length: " << TAGW << std::endl;
    std::cout << "TAGE history lengths: ";
    for (u64 i=0; i<NUMG; i++) std::cout << gfolds.HLEN[i] << " ";
    std::cout << std::endl;
  }


  pred_output predict(val<64> inst_pc)
  {
    // executed once per cycle
    static bool first_cycle = true;

    // do buffered updates
    if (! first_cycle) {
      execute_if(val<1>{update_rank!=0}, [&](){
	update_cycle();
	update_cycle_history();
      });
    }
    first_cycle = false;
    update_rank = 0;

    lineaddr = inst_pc.fo1() >> loglinebytes;
    lineaddr.fanout(hard<1+NUMG*2>{});
    gfolds.fanout(hard<2>{});

    // compute indexes and hashed tags
    bindex = lineaddr;
    bindex.fanout(hard<lineinst>{});
    for (u64 i=0; i<NUMG; i++) {
      gindex[i] = lineaddr ^ gfolds.template get<0>(i);
      htag[i] = val<HTAGBITS>{lineaddr}.reverse() ^ gfolds.template get<1>(i);
    }
    gindex.fanout(hard<3>{});  
    htag.fanout(hard<2>{});

    // read tables
    for (u64 offset=0; offset<lineinst; offset++) {
      readb[offset] = bim[offset].read(bindex);
    }
    readb.fanout(hard<2>{});
    for (u64 i=0; i<NUMG; i++) {
      readt[i] = gtag[i].read(gindex[i]);
      readc[i] = gctr[i].read(gindex[i]);
      readu[i] = ubit[i].read(gindex[i]);
    }
    readt.fanout(hard<lineinst+1>{});
    readc.fanout(hard<3>{});
    readu.fanout(hard<2>{});
    notumask = ~readu.concat();
    notumask.fanout(hard<2>{});

    // gather prediction bits for each offset
    arr<val<1>,NUMG> gpreds_split = [&](int i) {return val<1>{readc[i]>>(CTR-1)};};
    val<NUMG> gpreds = gpreds_split.fo1().concat();
    gpreds.fanout(hard<lineinst>{});
    arr<val<NUMG+1>,lineinst> preds = [&](u64 offset){return concat(val<1>{readb[offset]>>1},gpreds);};
    preds.fanout(hard<2*lineinst>{});

    // hashed tags comparisons
    arr<val<1>,NUMG> htagcmp_split = [&](int i){return val<HTAGBITS>{readt[i]} == htag[i];};
    val<NUMG> htagcmp = htagcmp_split.fo1().concat();
    htagcmp.fanout(hard<lineinst>{});

    // generate match mask for each offset
    static_loop<lineinst>([&]<u64 offset>(){
      arr<val<1>,NUMG> tagcmp = [&](int i){return val<loglineinst>{readt[i]>>HTAGBITS} == hard<offset>{};};
      match[offset] = concat(val<1>{1}, tagcmp.fo1().concat() & htagcmp); // bimodal is default when no match
    });
    match.fanout(hard<2>{});

    // for each offset, find longest match and select primary prediction
    for (u64 offset=0; offset<lineinst; offset++) {
      match1[offset] = match[offset].one_hot();
    }
    match1.fanout(hard<3>{});
    for (u64 offset=0; offset<lineinst; offset++) {
      pred1[offset] = (match1[offset] & preds[offset]) != hard<0>{};
    }
    pred1.fanout(hard<2>{});

    // for each offset, find second longest match and select secondary prediction
    for (u64 offset=0; offset<lineinst; offset++) {
      match2[offset] = (match[offset]^match1[offset]).one_hot();
    }
    match2.fanout(hard<2>{});
    for (u64 offset=0; offset<lineinst; offset++) {
      pred2[offset] = (match2[offset] & preds[offset]) != hard<0>{};
    }
    pred2.fanout(hard<2>{});

#ifdef USE_META
    meta.fanout(hard<2>{});
    arr<val<1>,NUMG> weakctr = [&](int i) {return (readc[i]==hard<WEAK0>{}) | (readc[i]==hard<WEAK1>{});};
    val<NUMG> coldctr = notumask & weakctr.fo1().concat();
    coldctr.fanout(hard<lineinst>{});
    val<1> metasign = (meta[METAPIPE-1] >= hard<0>{});
    metasign.fanout(hard<lineinst>{});
    for (u64 offset=0; offset<lineinst; offset++) {
      newly_alloc[offset] = (match1[offset] & coldctr) != hard<0>{};
      newly_alloc.fanout(hard<2>{});
      arr<val<1>,3> altselconds = {metasign, newly_alloc[offset], match2[offset]!=hard<0>{}};
      prediction[offset] = select(altselconds.fo1().fold_and(),pred2[offset],pred1[offset]);
    }
#else
    prediction = pred1;
#endif
    prediction.fanout(hard<2>{});
    arr<val<1>,lineinst> valid = [](){return 1;};
    return {prediction,valid.fo1()};
  }


  void update_cycle()
  {
    // executed once per cycle
    prediction.fanout(hard<2>{});
    bindex.fanout(hard<lineinst>{});
    gindex.fanout(hard<3>{});
    htag.fanout(hard<3>{});
    readb.fanout(hard<2>{});
    readt.fanout(hard<4>{});
    readc.fanout(hard<2>{});
    notumask.fanout(hard<2>{});
    match1.fanout(hard<2>{});
    match2.fanout(hard<2>{});
    pred1.fanout(hard<2>{});
    pred2.fanout(hard<2+NUMG>{});
    update_offset.fanout(hard<lineinst+NUMG+1>{});
    update_dir.fanout(hard<3>{});
    update_nextinst.fanout(hard<2>{});
    lineaddr.fanout(hard<2>{});
    gfolds.fanout(hard<2>{});
#ifdef USE_META
    meta.fanout(hard<2>{});
    newly_alloc.fanout(hard<2>{});
#endif
    u64 update_valid = (u64(1)<<update_rank)-1;

    if (update_rank==0) {
      last_offset = lineinst-1;
      last_dir = 0;
    } else {
      last_offset = update_offset[update_rank-1];
      last_dir = update_dir[update_rank-1];
    }
    last_dir.fanout(hard<2>{});
    last_offset.fanout(hard<lineinst+4*NUMG>{});

    static_loop<lineinst>([&]<u64 offset>{
      arr<val<1>,lineinst> match_offset = [&](u64 i){return update_offset[i] == hard<offset>{};};
      update_mask[offset] = match_offset.fo1().concat() & update_valid;
    });
    update_mask.fanout(hard<2>{});

    arr<val<1>,lineinst> is_branch = [&](u64 offset){
      return update_mask[offset] != hard<0>{};
    };
    is_branch.fanout(hard<3>{});

    val<lineinst> actualdirs = update_dir.concat();
    actualdirs.fanout(hard<lineinst+NUMG>{});

    arr<val<1>,lineinst> branch_dir = [&](u64 offset){
      return (actualdirs & update_mask[offset]) != hard<0>{};
    };
    branch_dir.fanout(hard<2>{});

    // at most one branch, the last one, can be mispredicted (as a misprediction ends a block)
    static_loop<lineinst>([&]<u64 offset>{
      update_last_branch[offset] = (last_offset == hard<offset>{});
    });
    update_last_branch.fanout(hard<2>{});
    val<1> last_pred = ((update_last_branch.concat() & prediction.concat()) != hard<0>{});
    val<1> mispred = (last_pred.fo1() != last_dir);
    mispred.fanout(hard<2*NUMG+1>{});

#ifdef USE_META
    // update meta counter
    arr<val<1>,lineinst> altdiff = [&](u64 offset){
      // for each offset, tell if primary and secondary predictions differ
      return (match2[offset] != hard<0>{}) & (pred2[offset] != pred1[offset]);
    };
    arr<val<2,i64>,lineinst> meta_incr = [&](u64 offset) -> val<2,i64> {
      val<1> update_meta = is_branch[offset] & altdiff[offset].fo1() & newly_alloc[offset];
      val<1> bad_pred2 = (pred2[offset] != branch_dir[offset]);
      return select(update_meta.fo1(),concat(bad_pred2.fo1(),val<1>{1}),val<2>{0});
    };
    for (u64 i=METAPIPE-1; i!=0; i--) {
      meta[i] = meta[i-1];
    }
    auto newmeta = meta[0] + meta_incr.fo1().fold_add();
    newmeta.fanout(hard<3>{});
    using meta_t = valt<decltype(meta[0])>;
    meta[0] = select(newmeta>meta_t::maxval, meta_t{meta_t::maxval}, select(newmeta<meta_t::minval, meta_t{meta_t::minval}, meta_t{newmeta}));
#endif

    // select some candidate entries for allocation
    val<NUMG> mispmask = mispred.replicate(hard<NUMG>{}).concat();
    arr<val<1>,NUMG> last_tagcmp = [&](int i){return readt[i] == concat(last_offset,htag[i]);};
    val<NUMG+1> last_match1 = last_tagcmp.fo1().append(1).concat().one_hot();
    last_match1.fanout(hard<2>{});
    val<NUMG> postmask = mispmask.fo1() & val<NUMG>(last_match1-1);
    postmask.fanout(hard<2>{});
    val<NUMG> candallocmask = postmask & notumask; // candidate post entries for allocation
    candallocmask.fanout(hard<2>{});

    // if multiple candidate entries, we select a single one, with some randomization
    val<NUMG> collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    val<NUMG> collamask1 = collamask.one_hot();
    collamask1.fanout(hard<3>{});
    val<NUMG> collamask2 = (collamask^collamask1).one_hot();
    val<NUMG> collamask12 = select(val<2>{std::rand()}==hard<0>{}, collamask2.fo1(), collamask1);
    val<NUMG> allocmask = collamask12.fo1().reverse();
    allocmask.fanout(hard<4>{});
    arr<val<1>,NUMG> allocate = allocmask.make_array(val<1>{});
    allocate.fanout(hard<3>{});

    // update the bimodal counters that provided a primary prediction
    for (u64 offset=0; offset<lineinst; offset++) {
      val<1> bim_primary = match1[offset] >> NUMG;
      execute_if(is_branch[offset] & bim_primary.fo1(), [&](){
	bim[offset].write(bindex, update_ctr(readb[offset],branch_dir[offset]));
      });
    }

    // associate a branch direction to each global table
    arr<val<1>,NUMG> bdir = [&](u64 i) {
      val<loglineinst> tag_offset = readt[i] >> HTAGBITS;
      val<loglineinst> offset = select(allocate[i],last_offset,tag_offset.fo1());
      offset.fanout(hard<lineinst>{});
      arr<val<1>,lineinst> match_offset = [&](u64 j){return update_offset[j] == offset;};
      return (match_offset.fo1().concat() & update_valid & actualdirs) != hard<0>{};
    };
    bdir.fanout(hard<2>{});
    
    // update a global counter if it provided a primary prediction or is in the allocated entry
    arr<val<NUMG+1>,lineinst> actual_match1 = [&](u64 offset){
      return select(is_branch[offset],match1[offset],val<NUMG+1>{0});
    };
    val<NUMG> primarymask = actual_match1.fo1().fold_or();
    primarymask.fanout(hard<2>{});
    execute_if(primarymask | allocmask, [&](u64 i){
      val<CTR> ctrinit = select(bdir[i],val<CTR>{WEAK1},val<CTR>{WEAK0});
      val<CTR> updated_ctr = select(allocate[i], ctrinit.fo1(), update_ctr(readc[i],bdir[i]));
      gctr[i].write(gindex[i],updated_ctr.fo1());
    });

    // overwrite the tag in the allocated entry
    execute_if(allocmask, [&](u64 i){gtag[i].write(gindex[i],concat(last_offset,htag[i]));});

    // what follows is for updating u bits

    // associate to each global table a bit telling if local prediction differs from secondary prediction
    arr<val<1>,NUMG> altdiffer = [&](u64 i){
      val<1> thispred = readc[i] >> (CTR-1);
      val<loglineinst> tag_offset = readt[i] >> HTAGBITS;
      return thispred.fo1() != pred2.select(tag_offset.fo1());
    };

    // associate to each global table a bit telling if prediction for this branch is correct
    arr<val<1>,NUMG> goodpred = [&](u64 i){
      val<loglineinst> tag_offset = readt[i] >> HTAGBITS;
      return select(tag_offset.fo1()==last_offset, ~mispred, val<1>{1});
    };

    // update the u bits
    val<NUMG> umask = primarymask & altdiffer.fo1().concat(); // u bits update mask
    // if all post entries have the u bit set, reset their u bits
    val<1> noalloc = (candallocmask == hard<0>{});
    val<NUMG> uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
    uclearmask.fanout(hard<2>{});
    arr<val<1>,NUMG> uclear = uclearmask.make_array(val<1>{});
    execute_if(umask.fo1() | allocmask | uclearmask, [&](u64 i) {
      val<1> newu = select(allocate[i] | uclear[i].fo1(), val<1>{0}, goodpred[i].fo1());
      ubit[i].write(gindex[i],newu.fo1());
    });

#ifdef RESET_UBITS
    uctr.fanout(hard<3>{});
    val<NUMG> allocmask1  = collamask1.reverse();
    allocmask1.fanout(hard<2>{});
    val<1> faralloc = (((last_match1>>3) | allocmask1).one_hot() ^ allocmask1) == hard<0>{};
    val<1> uctrsat = (uctr == hard<decltype(uctr)::maxval>{});
    uctrsat.fanout(hard<2>{});
    uctr = select(~mispred,uctr,select(uctrsat,val<decltype(uctr)::size>{0},update_ctr(uctr,faralloc.fo1())));
    execute_if(uctrsat,[&](){for (auto &uram : ubit) uram.reset();});
#endif
  }

  void update_cycle_history()
  {
    // update global history and folds with the address of the next block
    lineaddr.fanout(hard<2>{});
    last_dir.fanout(hard<2>{});
    if (update_rank == 0) {
      target_inst = 0;
    } else {
      target_inst = update_nextinst[update_rank-1];
    }
    target_inst.fanout(hard<2>{});
    val<PATHBITS> fall_through_inst = concat(val<PATHBITS-loglineinst>{lineaddr}+1,val<loglineinst>{0});
    val<PATHBITS> next_inst = select(~last_dir,fall_through_inst.fo1(),target_inst);
    gfolds.update(next_inst.fo1());
  }
  
  void update(val<64> branch_pc, val<1> dir, val<64> next_pc)
  {
    // on every conditional branch (update buffered)
    assert(update_rank < lineinst);
    update_offset[update_rank] = branch_pc.fo1() >> loginstbytes;
    update_dir[update_rank] = dir.fo1();
    update_nextinst[update_rank] = next_pc.fo1() >> loginstbytes;
    update_rank++;
  }

};
