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
  //synthetic_trace<> strace{14,0.001};
  synthetic_trace<1000,128,1024,5> strace{99,0.0};
  
  auto next_branch()
  {
    auto [branch_pc,branch_dir,nextpc] = strace.next();
    hcm::val<64> pc = branch_pc;
    hcm::val<1> dir = branch_dir;
    pc.set_time(t);
    dir.set_time(t); // FIXME?
    t += hcm::panel.clock_cycle_ps;
    return std::tuple{pc,dir};
  }

public:

  simulator()
  {
    hcm::panel.clock_cycle_ps = 250;
  }

  void run(predictor &p, int trace_length=10)
  {
    for (int i=0; i<trace_length; i++) {
      nbranch++;
      auto [pc,dir] = next_branch();
      hcm::val<1> pred = p.predict(pc);
      assert(pred.time() >= pc.time());
      uint64_t latency_ps = pred.time()-pc.time();
      if (latency_ps > max_pred_lat_ps)
	max_pred_lat_ps = latency_ps;
      if (pred.get() != dir.get()) {
	nmisp++;
      }
      p.update(pc,dir);
      hcm::panel.next_cycle();
    }
  }

  ~simulator()
  {
    strace.print_stats("synthetic trace ");
    std::cout << "branches: " << nbranch << std::endl;
    std::cout << "mispredicted: " << nmisp << std::endl;
    if (nbranch!=0)
      std::cout << "mispredict ratio: " << f64(nmisp)/nbranch << std::endl;
    hcm::panel.print();
    std::cout << "max prediction latency (cycle): " << double(max_pred_lat_ps) / hcm::panel.clock_cycle_ps << std::endl;
  }
};
