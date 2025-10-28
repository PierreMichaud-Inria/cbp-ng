#include "predictors/bimodal.hpp"
#include "predictors/gshare.hpp"
#include "predictors/tage.hpp"

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
//using branch_predictor = bimodal<17,6>;
//using branch_predictor = gshare<17,12,6>;
using branch_predictor = tage<8,11,12,11,100,6>;
#endif
