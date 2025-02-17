#include "gshare.hpp"
#include "tage.hpp"

//using branch_predictor = gshare<12>;
using branch_predictor = tage<8,12,14,8,200>;
