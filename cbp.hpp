#include <cstdint>
#include <iostream>
#include <tuple>
#include "harcom.hpp"

//using namespace hcm;


const std::tuple<uint64_t,bool> trace[] = {
  {0x1000,0},
  {0x1004,1},
  {0xFFFF0000,0},
  {0xFFFF0004,1},
  {0x1008,1}
};


struct predictor {
  virtual hcm::val<1> predict(hcm::val<64> pc) = 0;
  virtual void update(hcm::val<64> pc, hcm::val<1> dir) = 0;
  virtual ~predictor() {}
};


class simulator {
  uint64_t n = 0;
  uint64_t nbranch = 0;
  uint64_t nmisp = 0;
  uint64_t max_pred_lat_ps = 0;
  uint64_t t = 0; // ps

  template<uint64_t N>
  auto next_branch(const std::tuple<uint64_t,bool> (&trace)[N])
  {
    std::tuple<hcm::val<64>,hcm::val<1>> branch = trace[n];
    std::get<0>(branch).set_time(t);
    std::get<1>(branch).set_time(t); // FIXME?
    n++;
    t += hcm::panel.clock_cycle_ps;
    if (n==N) n = 0;
    return branch;
  }

public:

  simulator()
  {
    hcm::panel.clock_cycle_ps = 250;
    hcm::panel.clock_cycle_ps.print("clock cycle (ps): ");
    std::cout << std::setprecision(3);
    std::cout << "clock frequency (GHz): " << 1000./hcm::panel.clock_cycle_ps << std::endl;
  }
  
  void run(predictor &p, int trace_length=10)
  {
    for (int i=0; i<trace_length; i++) {
      nbranch++;
      auto [pc,dir] = next_branch(trace);
      hcm::val<1> pred = p.predict(pc);
      assert(pred.time() >= pc.time());
      uint64_t latency_ps = pred.time()-pc.time();
      if (latency_ps > max_pred_lat_ps)
	max_pred_lat_ps = latency_ps;
      if (pred.get() != dir.get()) {
	nmisp++;
      }
      p.update(pc,dir);
    }
  }

  ~simulator()
  {
    std::cout << "branches: " << nbranch << std::endl;
    std::cout << "mispredicted: " << nmisp << std::endl;
    std::cout << std::setprecision(3);
    hcm::panel.storage.print("storage (bits): ");
    std::cout << "max prediction latency (cycle): " << double(max_pred_lat_ps) / hcm::panel.clock_cycle_ps << std::endl;
    if (t!=0) {
      std::cout << "power (mW): " << hcm::panel.energy_fJ / t << std::endl;
    }
  }
};
