#include "cbp.hpp"
#include "branch_predictor.hpp"


simulator sim;
branch_predictor pred;


int main()
{
  sim.run(pred,1000000);
}
