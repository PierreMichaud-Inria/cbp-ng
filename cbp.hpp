#include <cstdint>
#include <iostream>
#include <tuple>
#include "synthetic_trace.hpp"
#include "harcom.hpp"

using namespace hcm;


static constexpr u64 loglinebytes = 6;
static constexpr u64 linebytes = 1 << loglinebytes; 
static constexpr u64 loginstbytes = 2;
static constexpr u64 instbytes = 1 << loginstbytes;
static constexpr u64 loglineinst = loglinebytes - loginstbytes;
static constexpr u64 lineinst = 1 << loglineinst;


u64 inst_offset(u64 inst_pc)
{
  return (inst_pc >> loginstbytes) % lineinst;
}


struct predictor {
  using predbits = arr<val<1>,lineinst>;
  using validbits = arr<val<1>,lineinst>;
  using pred_output = std::tuple<predbits,validbits>;
  virtual pred_output predict(val<64> inst_pc) = 0;
  virtual void update(val<64> branch_pc, val<1> dir, val<64> next_pc) = 0;
  virtual void update_cycle() = 0;
};


class harcom_superuser {
  predictor &pred;
  u64 ninst = 0;
  u64 ncycle = 0;
  u64 nbranch = 0;
  u64 nmisp = 0;
  u64 max_pred_lat_ps = 0;
  u64 t = 0; // ps

  auto next_instruction(synthetic_trace &trace)
  {
    auto [inst_pc, is_branch, is_taken, next_pc] = trace.next_instruction();
    ninst++;
    if (is_branch) nbranch++;
    return std::tuple {inst_pc, is_branch, is_taken, next_pc};
  }
  
public:

  harcom_superuser(predictor &pred) : pred(pred)
  {
    panel.clock_cycle_ps = 250;
    panel.make_floorplan();
  }

  void run(synthetic_trace &trace, u64 trace_length=10)
  {
    trace_length += nbranch;
    bool pred_available = false;
    std::array<u64,lineinst> taken;
    std::array<u64,lineinst> valid;

    while (nbranch < trace_length) {

      // one instruction per iteration
      auto [pc, is_branch, is_taken, next_pc] = next_instruction(trace);

      if (! pred_available) {
	// get predictions for new block
	auto [preds,valids] = pred.predict({pc,t});
	taken = preds.fo1().get();
	valid = valids.fo1().get();
	pred_available = true;
	u64 predtime = 0;
	for (u64 i=0; i<lineinst; i++) {
	  predtime = std::max(predtime,preds[i].time());
	  predtime = std::max(predtime,valids[i].time());
	}
	if (predtime > t) {
	  u64 latency_ps = predtime - t;
	  if (latency_ps > max_pred_lat_ps)
	    max_pred_lat_ps = latency_ps;
	}
      }

      assert(valid[inst_offset(pc)]);
      bool mispredict = is_branch && (is_taken != taken[inst_offset(pc)]);
      if (mispredict) {
	nmisp++;
      }
      if (is_branch) {
	// immediate update (not realistic)
	pred.update({pc,t},{is_taken,t},{next_pc,t});
      }

      bool is_jump = (next_pc != (pc + instbytes));
      bool line_end = (inst_offset(pc) == (lineinst-1));
      bool last_pred = line_end || ! valid[inst_offset(pc+instbytes)];
      bool end_of_trace = (nbranch == trace_length);
      if (is_jump || line_end || mispredict || last_pred || end_of_trace) {
	// end of block
	pred.update_cycle();
	panel.next_cycle();
	t += hcm::panel.clock_cycle_ps;
	ncycle++;
	pred_available = false; // need new prediction
      }
    }
  }

  ~harcom_superuser()
  {
    std::cout << "instructions: " << ninst << std::endl;
    std::cout << "cycles: " << ncycle << std::endl;
    std::cout << "branches: " << nbranch << std::endl;
    std::cout << "mispredicted: " << nmisp << std::endl;
    if (nbranch!=0)
      std::cout << "mispredict ratio: " << f64(nmisp)/nbranch << std::endl;
    panel.print();
    std::cout << "max prediction latency (cycle): " << double(max_pred_lat_ps) / panel.clock_cycle_ps << std::endl;
  }
};

