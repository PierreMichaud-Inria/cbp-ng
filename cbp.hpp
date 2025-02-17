#include <tuple>
#include "harcom.hpp"

using namespace hcm;

inline constexpr u64 clock_cycle_ps = 250;


const std::tuple<u64,bool> trace[] = {
  {0x1000,0},
  {0x1004,1},
  {0xFFFF0000,0},
  {0xFFFF0004,1},
  {0x1008,1}
};


struct predictor {
  virtual val<1> predict(val<64> pc) = 0;
  virtual void update(val<64> pc, val<1> dir) = 0;
  virtual ~predictor() {}
};


class simulator {
  u64 n = 0;
  u64 nbranch = 0;
  u64 nmisp = 0;
  u64 max_pred_lat_ps = 0;
  u64 t = 0; // ps

  template<u64 N>
  auto next_branch(const std::tuple<u64,bool> (&trace)[N])
  {
    std::tuple<val<64>,val<1>> branch = trace[n];
    std::get<0>(branch).set_time(t);
    std::get<1>(branch).set_time(t); // FIXME?
    n++;
    t += clock_cycle_ps;
    if (n==N) n = 0;
    return branch;
  }

public:
  
  void run(predictor &p, int trace_length=10)
  {
    for (int i=0; i<trace_length; i++) {
      nbranch++;
      auto [pc,dir] = next_branch(trace);
      val<1> pred = p.predict(pc);
      assert(pred.time() >= pc.time());
      u64 latency_ps = pred.time()-pc.time();
      max_pred_lat_ps = max(max_pred_lat_ps,latency_ps);
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
    std::cout << std::setprecision(4);
    metrics.storage.print("storage (bits): ");
    std::cout << "max prediction latency (cycle): " << f64(max_pred_lat_ps) / clock_cycle_ps << std::endl;
    if (t!=0) {
      std::cout << "power (mW): " << metrics.energy_fJ / t << std::endl;
    }
  }
};
