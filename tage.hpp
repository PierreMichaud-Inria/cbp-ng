#include <cstdlib>
#include "harcom.hpp"

using namespace hcm;


template<u64 MAXL>
struct path_history {
  static constexpr u64 N = MAXL / 64;
  static constexpr u64 R = MAXL % 64;

  arr<reg<64>,N> h;
  reg<R> hr;
  
  void update(ival auto in)
  {
    if constexpr (N==0) {
      hr = (hr << 1) ^ in;
    } else {
      if constexpr (R!=0) {
	hr = concat(val<R-1>(hr),h[N-1][63]);
      }
      for (int i=N-1; i>=1; i--) {
	h[i] = concat(val<63>(h[i]),h[i-1][63]);
      }
      h[0] = (h[0] << 1) ^ in;
    }
  }

  template<u64 STEP, u64 LENGTH>
  auto make_array()
  {
    if constexpr (R==0) {
      return h.template make_array<STEP,LENGTH>();
    } else {
      return h.append(hr).template make_array<STEP,LENGTH>();
    }
  }
};


template<u64 N>
auto update_ctr(val<N> ctr, val<1> incr)
{
  val<N> ctrmax = -1;
  val<N> incsat = select(ctr==ctrmax,ctr,ctr+1);
  val<N> decsat = select(ctr==0,ctr,ctr-1);
  return select(incr,incsat,decsat);
}


template<u64 NUMG, u64 LOGG, u64 LOGB, u64 TAGW, u64 GHIST>
struct tage : predictor {
  static constexpr u64 CTR = 3;
  static constexpr u64 MINHIST = 3;
  static constexpr u64 NG = 1<<LOGG;
  static constexpr u64 NB = 1<<LOGB;
  
  ram<val<2>,NB> bim; // bimodal table
  ram<val<TAGW>,NG> gtag[NUMG]; // global tables tags
  ram<val<CTR>,NG> gctr[NUMG]; // global tables counters
  ram<val<1>,NG> ubit[NUMG]; // "useful" bits
  path_history<GHIST> ph;
  reg<LOGB> bi; // bimodal table index
  arr<reg<LOGG>,NUMG> gi; // global tables indexes
  arr<reg<TAGW>,NUMG> gt; // computed tags
  reg<2> readb; // read bimodal counter
  arr<reg<TAGW>,NUMG> readt;  // read tags
  arr<reg<CTR>,NUMG> readc; // read global counters
  arr<reg<1>,NUMG> readu; // read u bits
  arr<reg<1>,NUMG> tagcmp; // results of tag comparisons
  reg<NUMG+1> match;
  reg<NUMG+1> match1; // longest match
  reg<NUMG+1> match2; // second longest match
  reg<1> prediction;

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
    // TODO: use altpred when low conf prediction
    // compute indexes
    bi = pc;
    static_loop<NUMG> ([&]<int I>() {
	auto a = ph.template make_array<LOGG,HLEN[I]>();
	a.append(pc);
	gi[I] = a.xor_all();
    });
    // compute tags
    static_loop<NUMG> ([&]<int I>() {
	auto a = ph.template make_array<TAGW,HLEN[I]>();
	a.append(val<TAGW>(pc).reverse());
	gt[I] = a.xor_all();
    });
    // read tables
    readb = bim.read(bi);
    for (int i=0; i<NUMG; i++) {
      readt[i] = gtag[i].read(gi[i]);
      readc[i] = gctr[i].read(gi[i]);
      readu[i] = ubit[i].read(gi[i]);
    }
    // extract prediction bits
    arr<val<1>,NUMG> gpreds = [&](int i) {return readc[i][CTR-1];};
    val<NUMG+1> preds = gpreds.append(readb[1]).concat();
    // compare tags
    for (int i=0; i<NUMG; i++) {
      tagcmp[i] = (gt[i] == readt[i]);
    }
    // find the longest match
    match = tagcmp.append(1).concat();
    match1 = match.priority_encode();
    // final prediction
    prediction = (match1 & preds) != 0;
    return prediction;
  }

  void update(val<64> pc, val<1> dir)
  {
    // TODO: periodic reset of u bits
    val<1> goodpred = (prediction == dir);
    val<1> mispred = (prediction != dir);
    // update the counter that provided the prediction
    auto update_counter = [&] (int i) {
      if (i<NUMG) {
	gctr[i].write(gi[i],update_ctr(readc[i],dir));
      } else {
	bim.write(bi,update_ctr(readb,dir));
      }
    };
    execute(match1,update_counter);
    // determine the alternate prediction
    match2 = (match^match1).priority_encode(); // 2nd longest match
    arr<val<1>,NUMG> gpreds = [&](int i) {return readc[i][CTR-1];};
    val<NUMG+1> preds = gpreds.append(readb[1]).concat();
    val<1> altpred = (match2 & preds) != 0;
    val<1> altdiff = (altpred != prediction);
    // update u bit of providing entry
    val<NUMG> umask = match1 & altdiff.make_array<NUMG>().concat();
    execute(umask,[&](int i){ubit[i].write(gi[i],goodpred);});
    val<NUMG> postmask = match1-1; // histories longer than providing entry ("post" entries)
    val<NUMG> candallocmask = postmask & ~readu.concat(); // candidate post entries for allocation
    // if multiple candidate entries, we select a single one
    val<NUMG> collamask = candallocmask.reverse();
    val<NUMG> collamask1 = collamask.priority_encode();
    val<NUMG> collamask2 = (collamask^collamask1).priority_encode();
    val<2> rand2 = rand();
    val<NUMG> collamask12 = select(rand2==0,collamask2,collamask1);
    val<NUMG> allocmask = collamask12.reverse();
    auto allocate = [&] (int i) {
      gtag[i].write(gi[i],gt[i]);
      ubit[i].write(gi[i],0);
      val<CTR> initctr = dir.make_array<CTR>().concat() ^ ((1<<(CTR-1))-1);
      gctr[i].write(gi[i],initctr);
    };
    execute(allocmask,allocate);
    // if all post entries have the u bit set, reset their u bits
    val<1> noalloc = (candallocmask == 0);
    val<NUMG> uclearmask = postmask & noalloc.make_array<NUMG>().concat();
    execute(uclearmask,[&](int i){ubit[i].write(gi[i],0);});
    // update path history
    ph.update(concat(val<5>(pc),dir));
  }
  
};
