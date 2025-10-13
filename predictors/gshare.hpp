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
  arr<reg<loglineinst>,lineinst> update_offset;
  arr<reg<1>,lineinst> update_dir;
  arr<reg<GHIST>,lineinst> update_nextinst;
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
    // do the PHT updates
    for (u64 offset=0; offset<lineinst; offset++) {
      auto [is_branch,dir] = get_info(offset);
      execute_if(is_branch,[&](){
	pht[offset].write(index,update_ctr(ctr[offset],dir));
      });
    }
    // update the global history with the address of the next block
    val<1> fall_through = (update_rank==0)? val<1>{1} : ~update_dir[update_rank-1];
    val<GHIST> fall_through_inst = concat(lineaddr+1,val<loglineinst>{0});
    val<GHIST> target_inst = (update_rank==0)? val<GHIST>{0} : update_nextinst[update_rank-1];
    val<GHIST> next_inst = select(fall_through,fall_through_inst,target_inst);
    ghist = (ghist << 1) ^ next_inst;
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
    lineaddr = inst_pc.fo1() >> loglinebytes;
    if constexpr (GHIST <= indexbits) {
      index = lineaddr ^ (ghist << (indexbits-GHIST));
    } else {
      index = ghist.make_array(val<indexbits>{}).append(lineaddr).fold_xor();
    }
    index.fanout(hard<lineinst*2>{});
    for (u64 i=0; i<lineinst; i++)
      ctr[i] = pht[i].read(index);
    static_assert(ctr.size == predbits::size);
    predbits pred = [&](u64 i){return val<1>(ctr[i]>>1);};
    validbits valid = [](){return 1;};
    return {pred.fo1(),valid.fo1()};
  }

  void update(val<64> branch_pc, val<1> dir, val<64> next_pc)
  {
    // on every conditional branch
    // update will be done in the next cycle, along with other updates in same line
    assert(update_rank < lineinst);
    update_offset[update_rank] = branch_pc >> loginstbytes;
    update_dir[update_rank] = dir;
    update_nextinst[update_rank] = next_pc >> loginstbytes;
    update_rank++;
  }
};
