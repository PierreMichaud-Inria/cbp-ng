#include "../harcom.hpp"
#include "common.hpp"

using namespace hcm;

// GEHL perceptron (see Seznec's ISCA 2005 paper)


#define DYNTHRES


template<u64 NUMG, u64 LOGG, u64 LOGB, u64 CTR, u64 GHIST>
struct gehl : predictor {
  static constexpr u64 NUMW = NUMG+1;
  static constexpr u64 MINHIST = 3;
  static constexpr u64 NG = 1<<LOGG;
  static constexpr u64 NB = 1<<LOGB;
  static constexpr u64 TCBITS = 7;
  static constexpr u64 OFFSET = NUMW/2;
  static constexpr u64 SUMBITS = CTR + std::bit_width(NUMW);
  static constexpr u64 THETABITS = 5;
  static_assert((OFFSET>>CTR)==0);

  geometric_folds<NUMG,MINHIST,GHIST,LOGG> gfolds;
  ram<val<CTR,int>,NB> bias; // bias weight table
  ram<val<CTR,int>,NG> weight[NUMG]; // other weight tables
  reg<THETABITS> theta = NUMW; // threshold
  reg<LOGB> bi; // bias weight table index
  arr<reg<LOGG>,NUMG> wi; // weight tables indexes
  arr<reg<CTR,int>,NUMW> readw; // read weights
  reg<1> prediction; // final prediction
  reg<SUMBITS,int> sum;
#ifdef DYNTHRES
  reg<7,int> tc; // threshold counter (see Seznec, ISCA 2005)
#endif

  gehl()
  {
    std::cout << "GEHL bias weights: " << NB << std::endl;
    std::cout << "GEHL global weights: " << NUMG << " x " << NG << std::endl;
    std::cout << "GEHL bits per weight: " << CTR << std::endl; 
    std::cout << "GEHL history lengths: ";
    for (u64 i=0; i<NUMG; i++) std::cout << gfolds.HLEN[i] << " ";
    std::cout << std::endl;
  }

  val<1> predict(val<64> pc)
  {
    pc.fanout(hard<NUMG+1>{});
    gfolds.fanout(hard<2>{});
    // compute indexes
    bi = pc>>2;
    bi.fanout(hard<2>{});
    for (u64 i=0; i<NUMG; i++) {
      wi[i] = pc>>2 ^ gfolds.get(i).reverse();
    }
    wi.fanout(hard<2>{});
    // read tables
    for (u64 i=0; i<NUMG; i++) {
      readw[i] = weight[i].read(wi[i]);
    }
    readw[NUMG] = bias.read(bi);
    readw.fanout(hard<2>{});
    // compute prediction
    sum = readw.append(OFFSET).fold_add();
    sum.fanout(hard<2>{});
    prediction = (sum >= hard<0>{});
    prediction.fanout(hard<2>{});
    return prediction;
  }

  void update([[maybe_unused]] val<64> pc, val<1> dir)
  {
    dir.fanout(hard<NUMG+3>{});
    theta.fanout(hard<3>{});
    val<1> mispred = (prediction != dir);
    mispred.fanout(hard<4>{});
    val<1> weaksum = (absolute_value(sum) <= val<SUMBITS>{theta});
    weaksum.fanout(hard<2>{});
    val<1> upd = mispred | weaksum;
    upd.fanout(hard<NUMG+1>{});
    // update weights
    execute_if(upd,[&](){bias.write(bi,update_ctr(readw[NUMG],dir));});
    auto updmask = upd.replicate(hard<NUMG>{}).concat();
    execute_if(updmask.fo1(),[&](u64 i){weight[i].write(wi[i],update_ctr(readw[i],dir));});
    // update global history and folds
    auto branchbits = concat(val<5>{pc.fo1()>>2},dir);
    gfolds.update(branchbits);
#ifdef DYNTHRES
    // update threshold (see Seznec, ISCA 2005)
    tc.fanout(hard<4>{});
    val<1> updtc = mispred | weaksum;
    updtc.fanout(hard<2>{});
    val<1> tcsat = tc==tc.maxval | tc==tc.minval;
    tcsat.fanout(hard<2>{});
    tc = select(updtc,select(tcsat,val<TCBITS>{0},update_ctr(tc,mispred)),tc);
    theta = select(updtc & tcsat,update_ctr(theta,mispred),theta);
#endif
  }
};


