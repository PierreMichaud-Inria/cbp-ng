#include "predictors/bimodal.hpp"
#include "predictors/gshare.hpp"
#include "predictors/tage.hpp"

#ifdef PREDICTOR
using branch_predictor = PREDICTOR;
#else
//using branch_predictor = bimodal<14>;
//using branch_predictor = gshare<18,12>;
using branch_predictor = tage<8,11,12,11,200>;
#endif
