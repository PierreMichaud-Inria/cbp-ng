#include <cstdlib>
#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


#define USE_META


template<u64 NUMG, u64 LOGG, u64 LOGB, u64 TAGW, u64 GHIST>
struct tage : predictor {
  static constexpr u64 CTR = 3;
  static constexpr u64 MINHIST = 3;
  static constexpr u64 NG = 1<<LOGG;
  static constexpr u64 NB = 1<<LOGB;
  static constexpr u64 WEAK1 = 1<<(CTR-1);
  static constexpr u64 WEAK0 = WEAK1-1;
  static constexpr u64 METABITS = 4;

  ram<val<2>,NB> bim; // bimodal table
  ram<val<TAGW>,NG> gtag[NUMG]; // global tables tags
  ram<val<CTR>,NG> gctr[NUMG]; // global tables counters
  ram<val<1>,NG> ubit[NUMG]; // "useful" bits

  geometric_folds<NUMG,MINHIST,GHIST,LOGG,TAGW> gfolds;
  
  reg<LOGB> bi; // bimodal table index 
  arr<reg<LOGG>,NUMG> gi; // global tables indexes
  arr<reg<TAGW>,NUMG> gt; // computed tags

  reg<2> readb; // read bimodal counter
  arr<reg<CTR>,NUMG> readc; // read global counters
  arr<reg<1>,NUMG> readu; // read u bits

  reg<NUMG+1> match1; // longest match
  reg<NUMG+1> match2; // second longest match

  reg<1> pred1; // primary prediction
  reg<1> pred2; // alternate prediction
  reg<1> prediction; // final prediction

#ifdef USE_META
  reg<METABITS> meta; // select between pred1 and pred2
  reg<1> newly_alloc;
#endif

  static constexpr auto HLEN = [] () {
    // geometric history lengths
    // table 0 is the rightmost table and has the longest history    
    std::array<u64,NUMG> hlen;
    u64 prevhl = 0;
    for (u64 i=0; i<NUMG; i++) {
      u64 hl = MINHIST * pow(f64(GHIST)/MINHIST,f64(i)/(NUMG-1));
      hl = std::max(prevhl+1,hl);
      hlen[NUMG-1-i] = hl;
      prevhl = hl;
    }
    return hlen;
  }();

  tage()
  {
    std::cout << "TAGE global entries: " << NUMG << " x " << NG << std::endl;
    std::cout << "TAGE bimodal entries: " << NB << std::endl;
    std::cout << "TAGE global entry tag length: " << TAGW << std::endl;
    std::cout << "TAGE history lengths: ";
    for (u64 i=0; i<NUMG; i++) std::cout << HLEN[i] << " ";
    std::cout << std::endl;
  }

  val<1> predict(val<64> pc)
  {
    pc.fanout(hard<1+NUMG*2>{});
    // compute indexes
    bi = pc;
    bi.fanout(hard<2>{});
    for (u64 i=0; i<NUMG; i++) {
      gi[i] = pc ^ gfolds.template get<0>(i);
      gt[i] = val<TAGW>(pc).reverse() ^ gfolds.template get<1>(i);
    }
    gi.fanout(hard<6>{});  
    gt.fanout(hard<2>{});
    // read tables
    arr<val<TAGW>,NUMG> tags = [&](int i){return gtag[i].read(gi[i]);};
    readb = bim.read(bi);
    readb.fanout(hard<2>{});
    for (u64 i=0; i<NUMG; i++) {
      readc[i] = gctr[i].read(gi[i]);
      readu[i] = ubit[i].read(gi[i]);
    }
    readc.fanout(hard<4>{});
    readu.fanout(hard<2>{});
    
    // extract prediction bits
    arr<val<1>,NUMG> gpreds = [&](int i) {return val<1>{readc[i]>>(CTR-1)};};
    auto preds = gpreds.fo1().append(val<1>{readb>>1}).concat();
    preds.fanout(hard<2>{});
    // compare tags
    arr<val<1>,NUMG> tagcmp = [&](int i){return gt[i] == tags[i].fo1();};
    // find longest match
    auto match = tagcmp.fo1().append(1).concat(); // bimodal is default when no match
    match.fanout(hard<2>{});
    match1 = match.one_hot();
    match1.fanout(hard<8>{});
    pred1 = (match1 & preds) != hard<0>{};
    pred1.fanout(hard<3>{});
    // find second longest match
    match2 = (match^match1).one_hot();
    match2.fanout(hard<2>{});
    pred2 = (match2 & preds) != hard<0>{};
    pred2.fanout(hard<3>{});
#ifdef USE_META
    meta.fanout(hard<2>{});
    arr<val<1>,NUMG> weakctr = [&](int i) {return (readc[i]==hard<WEAK0>{}) | (readc[i]==hard<WEAK1>{});};
    newly_alloc = (val<NUMG>(match1) & weakctr.fo1().concat() & ~readu.concat()) != hard<0>{};
    newly_alloc.fanout(hard<2>{});
    prediction = select(newly_alloc & val<1>{meta>>(METABITS-1)},pred2,pred1);
#else
    prediction = pred1;
#endif
    prediction.fanout(hard<2>{});
    return prediction;
  }

  void update(val<64> pc, val<1> dir)
  {
    // TODO: periodic reset of u bits
    dir.fanout(hard<4+NUMG*2>{});
    auto goodpred = (prediction == dir);
    goodpred.fanout(hard<1+NUMG>{});
    auto mispred = ~goodpred;
    mispred.fanout(hard<NUMG>{});
    auto altdiff = (match2 != 0) & (pred2 != pred1);
    altdiff.fanout(hard<1+NUMG>{});
#ifdef USE_META
    execute_if(altdiff & newly_alloc, [&](){meta = update_ctr(meta,pred2==dir);});
#endif
    auto mispmask = mispred.replicate(hard<NUMG>{}).concat();
    auto postmask = mispmask.fo1() & val<NUMG>(match1-1);
    postmask.fanout(hard<2>{});
    auto candallocmask = postmask & ~readu.concat(); // candidate post entries for allocation
    candallocmask.fanout(hard<2>{});
    // if multiple candidate entries, we select a single one
    auto collamask = candallocmask.reverse();
    collamask.fanout(hard<2>{});
    auto collamask1 = collamask.one_hot();
    collamask1.fanout(hard<2>{});
    auto collamask2 = (collamask^collamask1).one_hot();
    auto collamask12 = select(val<2>{rand()}==hard<0>{}, collamask2.fo1(), collamask1);
    auto allocmask = collamask12.fo1().reverse();
    allocmask.fanout(hard<3>{});
    auto match1_split = match1.make_array(val<1>{});
    match1_split.fanout(hard<2>{});
    // update the bimodal counter if it provided the prediction
    execute_if(match1_split[NUMG], [&](){bim.write(bi,update_ctr(readb,dir));});
    // update the global counter if it provided the prediction or is in the allocated entry
    execute_if(val<NUMG>{match1} | allocmask, [&](u64 i) {
      auto newgctr = select(match1_split[i], update_ctr(readc[i],dir), select(dir,val<CTR>{WEAK1},val<CTR>{WEAK0}));
      gctr[i].write(gi[i],newgctr.fo1());
    });
    // update the tag in the allocated entry
    execute_if(allocmask, [&](u64 i){gtag[i].write(gi[i],gt[i]);});
    // update the ubit
    auto umask = (val<NUMG>{match1} & altdiff.replicate(hard<NUMG>{}).concat()); // u bit of providing entry
    umask.fanout(hard<2>{});
    auto umask_split = umask.make_array(val<1>{});
    // if all post entries have the u bit set, reset their u bits
    auto noalloc = (candallocmask == hard<0>{});
    auto uclearmask = postmask & noalloc.fo1().replicate(hard<NUMG>{}).concat();
    execute_if(umask | allocmask | uclearmask.fo1(), [&](u64 i) {
      ubit[i].write(gi[i],select(umask_split[i].fo1(),goodpred,val<1>{0}));
    });
    // update global history and folds
    auto branchbits = concat(val<5>{pc.fo1()},dir);
    gfolds.update(branchbits);
  }

};
