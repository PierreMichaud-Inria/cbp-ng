#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;


template<u64 LOGN, u64 GHIST>
struct gshare : predictor {
  static_assert(LOGN > loglineinst);
  static constexpr u64 indexbits = LOGN-loglineinst;

  reg<GHIST> ghist;
  reg<1> true_block;
  reg<indexbits> index;
  reg<std::max(indexbits,GHIST)> lineaddr;
  arr<reg<2>,lineinst> ctr;
  arr<reg<1>,lineinst> prediction; // final prediction for each offset
  
  // simulation artifacts (hardware cost may not be real)
  u64 num_branch = 0;
  arr<reg<loglineinst>,lineinst> branch_offset;
  arr<reg<1>,lineinst> branch_dir;
  arr<reg<GHIST>,lineinst> branch_nextinst;

  ram<val<2>,(1<<indexbits)> pht[lineinst];

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
    for (u64 i=0; i<lineinst; i++) {
      ctr[i] = pht[i].read(index);
    }
    ctr.fanout(hard<2>{});
    for (u64 i=0; i<lineinst; i++) {
      prediction[i] = ctr[i] >> 1;
    }
    prediction.fanout(hard<2>{});
    validbits valid = [](){return 1;};
    return {prediction,valid.fo1()};
  }

  void update(val<64> branch_pc, val<1> dir, val<64> next_pc)
  {
    // on every conditional branch
    // update will be done in the next cycle, along with other updates in same line
    assert(num_branch < lineinst);
    branch_offset[num_branch] = branch_pc.fo1() >> loginstbytes;
    branch_dir[num_branch] = dir.fo1();
    branch_nextinst[num_branch] = next_pc.fo1() >> loginstbytes;
    num_branch++;
  }

  void update_cycle([[maybe_unused]] val<1> mispredict)
  {
    // once per cycle
    if (num_branch == 0) {
      execute_if(~true_block,[&](){
	// previous block ended prematurely because of a mispredicted not-taken branch
	val<GHIST> fallthru_inst = concat(val<GHIST>{lineaddr}+1,val<loglineinst>{0});
	ghist = (ghist << 1) ^ fallthru_inst.fo1();
	true_block = 1;
      });
      return; // stop here
    }

    branch_offset.fanout(hard<lineinst>{});
    branch_dir.fanout(hard<4>{});
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
    branch_taken.fanout(hard<2>{});

    execute_if(is_branch.fo1().concat(), [&](u64 i){
      pht[i].write(index,update_ctr(ctr[i],branch_taken[i]));
    });

    val<loglineinst> last_offset = branch_offset[num_branch-1];
    last_offset.fanout(hard<2>{});
    true_block = ~mispredict.fo1() | branch_dir[num_branch-1] | (last_offset == hard<lineinst-1>{});

    // update the global history if this is a true block
    execute_if(true_block, [&](){
      val<GHIST> fallthru_inst = concat(val<GHIST>{lineaddr}+1,val<loglineinst>{0});
      val<GHIST> target_inst = branch_nextinst[num_branch-1];
      val<GHIST> next_inst = select(branch_dir[num_branch-1],target_inst.fo1(),fallthru_inst.fo1());
      ghist = (ghist << 1) ^ next_inst;
    });
    num_branch = 0; // done
  }

};
