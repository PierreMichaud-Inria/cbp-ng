#include "harcom.hpp"

using namespace hcm;


template<u64 LOGN>
struct gshare : predictor {
  static constexpr u64 N = 1<<LOGN;
  reg<LOGN> ghist;
  ram<val<2>,N> pht;
  reg<LOGN> index;
  reg<2> ctr;
  
  val<2> update_2bc(val<2> ctr, val<1> incr)
  {
    val<2> incsat = select(ctr==3,ctr,ctr+1);
    val<2> decsat = select(ctr==0,ctr,ctr-1);
    return select(incr,incsat,decsat);
  }

  val<1> predict(val<64> pc)
  {
    index = val<LOGN>(pc) ^ ghist;
    ctr = pht.read(index);
    return val<1>(ctr>>1);
  }

  void update(val<64> pc, val<1> dir)
  {
    pht.write(index, update_2bc(ctr,dir));
    ghist = concat(val<LOGN-1>(ghist),dir);
  }
};
