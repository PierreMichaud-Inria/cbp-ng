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
  arr<reg<loglineinst>,lineinst> branch_offset;
  arr<reg<1>,lineinst> branch_dir;
  u64 num_branch = 0; // no hardware cost (simulation artifact)

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
    branch_offset[num_branch] = branch_pc >> loginstbytes;
    branch_dir[num_branch] = dir;
    num_branch++;
  }

  val<2> get_branch_info(u64 offset)
  {
    // tells if the instruction at a given offset is a conditional branch
    // and gives its direction    
    arr<val<1>,lineinst> match_offset = [&](u64 i){
      return (i<num_branch) & (branch_offset[i]==offset);
    };
    arr<val<1>,lineinst> dir = [&](u64 i){
      return match_offset[i] & branch_dir[i];
    };
    val<1> is_branch = match_offset.fold_or();
    val<1> is_taken = dir.fold_or();
    return concat(is_branch.fo1(),is_taken.fo1());
  };

  void update_cycle()
  {
    // once per cycle
    branch_offset.fanout(hard<lineinst>{});
    branch_dir.fanout(hard<lineinst>{});
    arr<val<2>,lineinst> branch = [&](u64 offset){
      return get_branch_info(offset);
    };
    for (u64 i=0; i<lineinst; i++) {
      auto [is_branch,dir] = split<1,1>(branch[i].fo1());
      execute_if(is_branch.fo1(),[&](){
        bht[i].write(index,update_ctr(ctr[i],dir.fo1()));
      });
    }
    num_branch = 0; // done
  }
};
