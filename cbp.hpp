#include <cstdint>
#include <iostream>
#include <tuple>
#include "synthetic_trace.hpp"
#include "harcom.hpp"

//using namespace hcm;


struct predictor {
  virtual hcm::val<1> predict(hcm::val<64> pc) = 0;
  virtual void update(hcm::val<64> pc, hcm::val<1> dir) = 0;
  virtual ~predictor() {}
};


class simulator {
  uint64_t nbranch = 0;
  uint64_t nmisp = 0;
  uint64_t max_pred_lat_ps = 0;
  uint64_t t = 0; // ps
  
  auto next_branch(synthetic_trace &trace)
  {
    auto [branch_pc,branch_dir,nextpc] = trace.next();
    std::tuple branch {branch_pc,branch_dir,t};
    t += hcm::panel.clock_cycle_ps;
    nbranch++;
    return branch;
  }

public:

  simulator()
  {
    hcm::panel.clock_cycle_ps = 250;
  }

  void run(predictor &p, synthetic_trace &trace, int trace_length=10)
  {
    for (int i=0; i<trace_length; i++) {
      auto [branch_pc,branch_dir,start_time] = next_branch(trace);
      auto [prediction,pred_time] = p.predict({branch_pc,start_time}).get_vt();
      assert(pred_time >= start_time);
      uint64_t latency_ps = pred_time - start_time;
      if (latency_ps > max_pred_lat_ps)
	max_pred_lat_ps = latency_ps;
      if (prediction != branch_dir) {
	nmisp++;
      }
      p.update({branch_pc,start_time},{branch_dir,start_time/*FIXME?*/});
      hcm::panel.next_cycle();
    }
  }

  ~simulator()
  {
    std::cout << "branches: " << nbranch << std::endl;
    std::cout << "mispredicted: " << nmisp << std::endl;
    if (nbranch!=0)
      std::cout << "mispredict ratio: " << f64(nmisp)/nbranch << std::endl;
    hcm::panel.print();
    std::cout << "max prediction latency (cycle): " << double(max_pred_lat_ps) / hcm::panel.clock_cycle_ps << std::endl;
  }
};
