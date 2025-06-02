#include "predictors/gshare.hpp"
#include "predictors/tage.hpp"
#include "predictors/gehl.hpp"

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
//using branch_predictor = gshare<18,14>;
//using branch_predictor = gehl<8,12,12,4,200>;
using branch_predictor = tage<8,11,12,10,200>;
#endif
