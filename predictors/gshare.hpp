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
      index = val<LOGN>(pc) ^ (ghist<<(LOGN-GHIST));
    } else {
      index = ghist.make_array(val<LOGN>{}).append(pc).fold_xor();
    }
    ctr = pht.read(index);
    return val<1>(ctr>>1);
  }

  void update([[maybe_unused]] val<64> pc, val<1> dir)
  {
    pht.write(index, update_ctr(ctr,dir));
    ghist = concat(val<LOGN-1>(ghist),dir);
  }
};
