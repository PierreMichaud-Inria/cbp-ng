#include "gshare.hpp"
#include "tage.hpp"

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
//using branch_predictor = gshare<18>;
//using branch_predictor = tage<8,12,14,11,200>;
using branch_predictor = tage<8,10,12,11,200>;
#endif
