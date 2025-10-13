#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


template<u64 LOGN>
struct bimodal : predictor {
  static_assert(LOGN > loglineinst);
  static constexpr u64 indexbits = LOGN-loglineinst;
  ram<val<2>,(1<<indexbits)> bht[lineinst];
  reg<indexbits> index;
  arr<reg<2>,lineinst> ctr;
  arr<reg<loglineinst>,lineinst> update_offset;
  arr<reg<1>,lineinst> update_dir;
  u64 update_rank = 0; // no hardware cost (simulation artifact)

  void update_cycle()
  {
    // once per cycle
    update_offset.fanout(hard<lineinst>{});
    update_dir.fanout(hard<lineinst>{});
    // function object that tells if the instruction at a given offset is a branch and,
    // if it is one, gives its direction
    auto get_info = [&](u64 offset) {
      arr<val<1>,lineinst> match_offset = [&](u64 i){return (i<update_rank) & (update_offset[i]==offset);};
      arr<val<1>,lineinst> dir = [&](u64 i){return match_offset[i] & update_dir[i];};
      return std::tuple {match_offset.fold_or(), dir.fold_or()};
    };
    // do the updates
    for (u64 offset=0; offset<lineinst; offset++) {
      auto [is_branch,dir] = get_info(offset);
      execute_if(is_branch,[&](){
        bht[offset].write(index,update_ctr(ctr[offset],dir));
      });
    }
    update_rank = 0; // done
  }

  pred_output predict(val<64> inst_pc)
  {
    // once per cycle
    static bool first_cycle = true;
    if (! first_cycle) {
      update_cycle();
    }
    first_cycle = false;
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
    assert(update_rank < lineinst);
    update_offset[update_rank] = branch_pc >> loginstbytes;
    update_dir[update_rank] = dir;
    update_rank++;
  }
};
