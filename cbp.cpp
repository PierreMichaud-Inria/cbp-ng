#include <string>
#include <vector>
#include "synthetic_trace.hpp"
#include "harcom.hpp"
#include "cbp.hpp"
#include "branch_predictor.hpp"

harcom_superuser sim;
branch_predictor pred;
std::vector<synthetic_trace> traces;


int main(int argc, char * argv[])
{
  traces.emplace_back(1000,128,1024,5,99,0.0);
  traces.emplace_back(1000,8,101,5,14,0.001);
  
  u64 simlen = 1000000;
  u64 period = simlen / traces.size();
  if (argc>1) {
    simlen = std::stoul(argv[1]);
  }
  if (argc>2) {
    period = std::stoul(argv[2]);
  }  
  if (simlen < period) {
    period = simlen;
  }

  for (u64 n=0, t=0; n<simlen; n+=period, t=(t+1)%traces.size()) {
    sim.run(pred,traces[t],period);
  }
}
