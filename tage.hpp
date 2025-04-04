#include <cstdlib>
#include "harcom.hpp"

using namespace hcm;


template<u64 MAXL, u64 NUMG>
struct path_history {
  static constexpr u64 N = MAXL / 64;
  static constexpr u64 R = MAXL % 64;

  arr<reg<64>,N> h;
  reg<R> hr;

  path_history() {}

  template<u64 FO>
  void fanout(hard<FO> x)
  {
    h.fanout(hard<FO>{});
    hr.fanout(hard<FO>{});
  }

  void update(auto in)
  {
    constexpr u64 K = in.size;
    if constexpr (N==0) {
      hr = (hr<<1) ^ in;
    } else {
      if constexpr (R==0) {
	auto a = h.shift_left(val<1>{0});
	h[0] = a[0] ^ in;
	for (u64 i=1; i<N; i++) h[i] = a[i];
      } else {
	auto a = h.append(hr).shift_left(val<1>{0});
	static_assert(a.size==N+1);
	h[0] = a[0] ^ in;
	for (u64 i=1; i<N; i++) h[i] = a[i];
	hr = a[N];
      }
    }
  }

  template<u64 STEP, u64 LENGTH>
  auto make_array()
  {
    static_assert(LENGTH<=MAXL);
    static_assert(STEP!=0 && STEP<=64);
    constexpr u64 chunks = (LENGTH+STEP-1) / STEP;
    constexpr u64 pos = LENGTH % STEP;
    constexpr u64 mask = (u64(1)<<pos)-1;
    if constexpr (R==0) {
      auto a = h.make_array(val<STEP>{});
      arr<val<STEP>,chunks> chunked = [&](int i){return (mask!=0 && i==chunks-1)? a[i] & mask: a[i];};
      return chunked;
    } else {
      auto a = h.append(hr).make_array(val<STEP>{});
      arr<val<STEP>,chunks> chunked = [&](int i){return (mask!=0 && i==chunks-1)? a[i] & mask: a[i];};
      return chunked;
    }
  }
};


template<u64 N>
[[nodiscard]] val<N> update_ctr(val<N> ctr, val<1> incr)
{
  ctr.fanout(hard<6>{});
  val<N> incsat = select(ctr==ctr.maxval,ctr,ctr+1);
  val<N> decsat = select(ctr==ctr.minval,ctr,ctr-1);
  return select(incr,incsat,decsat);
}


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

  path_history<GHIST,NUMG> ph;

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
    int prevhl = 0;
    for (int i=0; i<NUMG; i++) {
      f64 hl = MINHIST * pow(f64(GHIST)/MINHIST,f64(i)/(NUMG-1));
      hl = max(prevhl+1,hl);
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
    for (int i=0; i<NUMG; i++) std::cout << HLEN[i] << " ";
    std::cout << std::endl;
  }

  val<1> predict(val<64> pc)
  {
    ph.fanout(hard<1+NUMG*2>{});
    pc.fanout(hard<1+NUMG*2>{});
    // compute indexes
    bi = pc;
    //bi.fanout(hard<2>{});
    static_loop<NUMG> ([&]<int I>() {
	gi[I] = ph.template make_array<LOGG,HLEN[I]>().append(pc).fold_xor();
    });
    gi.fanout(hard<6>{});
    // compute tags
    static_loop<NUMG> ([&]<int I>() {
	gt[I] = ph.template make_array<TAGW,HLEN[I]>().append(val<TAGW>(pc).reverse()).fold_xor();
    });
    //gt.fanout(hard<2>{});
    // read tables
    arr<val<TAGW>,NUMG> tags = [&](int i){return gtag[i].read(gi[i]);};
    readb = bim.read(bi);
    //readb.fanout(hard<2>{});
    for (int i=0; i<NUMG; i++) {
      readc[i] = gctr[i].read(gi[i]);
      readu[i] = ubit[i].read(gi[i]);
    }
    readc.fanout(hard<4>{});
    //readu.fanout(hard<2>{});

    // extract prediction bits
    arr<val<1>,NUMG> gpreds = [&](int i) {return val<1>{readc[i]>>(CTR-1)};};
    auto preds = gpreds.append(val<1>{readb>>1}).concat();
    //preds.fanout(hard<2>{});
    // compare tags
    arr<val<1>,NUMG> tagcmp = [&](int i){return gt[i] == tags[i];};
    // find longest match
    auto match = tagcmp.append(1).concat(); // bimodal is default when no match
    //match.fanout(hard<2>{});
    match1 = match.priority_encode();
    match1.fanout(hard<8>{});
    pred1 = (match1 & preds) != 0;
    //pred1.fanout(hard<3>{});
    // find second longest match
    match2 = (match^match1).priority_encode();
    //match2.fanout(hard<2>{});
    pred2 = (match2 & preds) != 0;
    //pred2.fanout(hard<3>{});
#ifdef USE_META
    //meta.fanout(hard<2>{});
    arr<val<1>,NUMG> weakctr = [&](int i) {return readc[i]==WEAK0 | readc[i]==WEAK1;};
    newly_alloc = (val<NUMG>(match1) & weakctr.concat() & ~readu.concat()) != 0;
    //newly_alloc.fanout(hard<2>{});
    prediction = select(newly_alloc & val<1>{meta>>(METABITS-1)},pred2,pred1);
#else
    prediction = pred1;
#endif
    //prediction.fanout(hard<2>{});
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
    execute(altdiff & newly_alloc, [&](){meta = update_ctr(meta,pred2==dir);});
#endif
    auto mispmask = mispred.replicate(hard<NUMG>{}).concat();
    auto postmask = mispmask & val<NUMG>(match1-1);
    //postmask.fanout(hard<2>{});
    auto candallocmask = postmask & ~readu.concat(); // candidate post entries for allocation
    //candallocmask.fanout(hard<2>{});
    // if multiple candidate entries, we select a single one
    auto collamask = candallocmask.reverse();
    //collamask.fanout(hard<2>{});
    auto collamask1 = collamask.priority_encode();
    //collamask1.fanout(hard<2>{});
    auto collamask2 = (collamask^collamask1).priority_encode();
    auto collamask12 = select(val<2>{rand()}==0, collamask2, collamask1);
    auto allocmask = collamask12.reverse();
    //allocmask.fanout(hard<3>{});
    auto match1_split = match1.make_array(val<1>{});
    // update the bimodal counter if it provided the prediction
    execute(match1_split[NUMG], [&](){bim.write(bi,update_ctr(readb,dir));});
    // update the global counter if it provided the prediction or is in the allocated entry
    execute(val<NUMG>{match1} | allocmask, [&](u64 i) {
      auto newgctr = select(match1_split[i], update_ctr(readc[i],dir), select(dir,val<CTR>{WEAK1},val<CTR>{WEAK0}));
      gctr[i].write(gi[i],newgctr);
    });
    // update the tag in the allocated entry
    execute(allocmask, [&](u64 i){gtag[i].write(gi[i],gt[i]);});
    // update the ubit
    auto umask = (val<NUMG>{match1} & altdiff.replicate(hard<NUMG>{}).concat()); // u bit of providing entry
    //umask.fanout(hard<2>{});
    auto umask_split = umask.make_array(val<1>{});
    // if all post entries have the u bit set, reset their u bits
    auto noalloc = (candallocmask == 0);
    noalloc.fanout(hard<NUMG>{});
    auto uclearmask = postmask & noalloc.replicate(hard<NUMG>{}).concat();
    execute(umask | allocmask | uclearmask, [&](u64 i) {
      ubit[i].write(gi[i],select(umask_split[i],goodpred,val<1>{0}));
    });
    // update path history
    ph.update(concat(val<5>{pc},dir));
  }
  
};
