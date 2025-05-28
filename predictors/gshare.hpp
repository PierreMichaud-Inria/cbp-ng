#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


template<u64 LOGN, u64 GHIST=10>
struct gshare : predictor {
  static constexpr u64 N = 1<<LOGN;
  reg<GHIST> ghist;
  ram<val<2>,N> pht;
  reg<LOGN> index;
  reg<2> ctr;
  
  val<1> predict(val<64> pc)
  {
    if constexpr (GHIST <= LOGN) {
      index = val<LOGN>(pc.fo1()>>2) ^ (ghist<<(LOGN-GHIST));
    } else {
      index = ghist.make_array(val<LOGN>{}).append(pc.fo1()>>2).fold_xor();
    }
    index.fanout(hard<2>{});
    ctr = pht.read(index);
    return val<1>(ctr>>1);
  }

  void update([[maybe_unused]] val<64> pc, val<1> dir)
  {
    dir.fanout(hard<2>{});
    pht.write(index, update_ctr(ctr,dir));
    ghist = concat(val<LOGN-1>(ghist),dir);
  }
};
