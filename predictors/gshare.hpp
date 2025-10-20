#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


template<u64 LOGN, u64 GHIST>
struct gshare : predictor {
  static_assert(LOGN > loglineinst);
  static constexpr u64 indexbits = LOGN-loglineinst;
  ram<val<2>,(1<<indexbits)> pht[lineinst];
  reg<GHIST> ghist;
  reg<indexbits> index;
  reg<std::max(indexbits,GHIST)> lineaddr;
  arr<reg<2>,lineinst> ctr;
  arr<reg<loglineinst>,lineinst> branch_offset;
  arr<reg<1>,lineinst> branch_dir;
  arr<reg<GHIST>,lineinst> branch_nextinst;
  u64 num_branch = 0; // no hardware cost (simulation artifact)

  pred_output predict(val<64> inst_pc)
  {
    // once per cycle
    lineaddr = inst_pc.fo1() >> loglinebytes;
    if constexpr (GHIST <= indexbits) {
      index = lineaddr ^ (ghist << (indexbits-GHIST));
    } else {
      index = ghist.make_array(val<indexbits>{}).append(lineaddr).fold_xor();
    }
    index.fanout(hard<lineinst*2>{});
    for (u64 i=0; i<lineinst; i++)
      ctr[i] = pht[i].read(index);
    predbits pred = [&](u64 i){return val<1>(ctr[i]>>1);};
    validbits valid = [](){return 1;};
    return {pred.fo1(),valid.fo1()};
  }

  void update(val<64> branch_pc, val<1> dir, val<64> next_pc)
  {
    // on every conditional branch
    // update will be done in the next cycle, along with other updates in same line
    assert(num_branch < lineinst);
    branch_offset[num_branch] = branch_pc >> loginstbytes;
    branch_dir[num_branch] = dir;
    branch_nextinst[num_branch] = next_pc >> loginstbytes;
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
	pht[i].write(index,update_ctr(ctr[i],dir));
      });
    }
    // update the global history with the address of the next block
    val<1> fall_through = (num_branch==0)? val<1>{1} : ~branch_dir[num_branch-1];
    val<GHIST> fall_through_inst = concat(lineaddr+1,val<loglineinst>{0});
    val<GHIST> target_inst = (num_branch==0)? val<GHIST>{0} : branch_nextinst[num_branch-1];
    val<GHIST> next_inst = select(fall_through,fall_through_inst,target_inst);
    ghist = (ghist << 1) ^ next_inst;
    num_branch = 0; // done
  }

};
