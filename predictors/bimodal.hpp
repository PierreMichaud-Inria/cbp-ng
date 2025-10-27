#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


template<u64 LOGN>
struct bimodal : predictor {
  static_assert(LOGN > loglineinst);
  static constexpr u64 indexbits = LOGN-loglineinst;

  reg<indexbits> index;
  arr<reg<2>,lineinst> ctr;

  // simulation artifacts (hardware cost may not be real)
  u64 num_branch = 0;
  arr<reg<loglineinst>,lineinst> branch_offset;
  arr<reg<1>,lineinst> branch_dir;

  ram<val<2>,(1<<indexbits)> bht[lineinst];

  pred_output predict(val<64> inst_pc)
  {
    // once per cycle
    index = inst_pc.fo1() >> loglinebytes;
    index.fanout(hard<lineinst*2>{});
    for (u64 offset=0; offset<lineinst; offset++)
      ctr[offset] = bht[offset].read(index);
    predbits pred = [&](u64 offset){return val<1>(ctr[offset]>>1);};
    validbits valid = [](){return 1;};
    return {pred.fo1(),valid.fo1()};
  }

  void update(val<64> branch_pc, val<1> dir, [[maybe_unused]] val<64> next_pc)
  {
    // on every conditional branch
    // update will be done in the next cycle, along with other updates in same line
    assert(num_branch < lineinst);
    branch_offset[num_branch] = branch_pc.fo1() >> loginstbytes;
    branch_dir[num_branch] = dir.fo1();
    num_branch++;
  }

  void update_cycle([[maybe_unused]] val<1> mispredict, [[maybe_unused]] val<64> next_pc)
  {
    // once per cycle
    if (num_branch == 0) {
      return;
    }
    branch_offset.fanout(hard<lineinst>{});
    branch_dir.fanout(hard<2>{});
    u64 update_valid = (u64(1)<<num_branch)-1;
    arr<val<lineinst>,lineinst> update_mask = [&](u64 offset){
      arr<val<1>,lineinst> match_offset = [&](u64 i){return branch_offset[i] == offset;};
      return match_offset.fo1().concat() & update_valid;
    };
    update_mask.fanout(hard<2>{});
    arr<val<1>,lineinst> is_branch = [&](u64 offset){
      return update_mask[offset] != hard<0>{};
    };
    val<lineinst> actualdirs = branch_dir.concat();
    actualdirs.fanout(hard<lineinst>{});
    arr<val<1>,lineinst> branch_taken = [&](u64 offset){
      return (actualdirs & update_mask[offset]) != hard<0>{};
    };
    // update predictor here
    execute_if(is_branch.fo1().concat(), [&](u64 i){
      bht[i].write(index,update_ctr(ctr[i],branch_taken[i].fo1()));
    });
    num_branch = 0; // done
  }
};
