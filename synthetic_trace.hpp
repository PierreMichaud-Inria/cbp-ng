#ifndef SYNTHETIC_TRACE_H
#define SYNTHETIC_TRACE_H

#include <cstdint>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>
#include <tuple>
#include <vector>
#include <array>
#include <bit>
#include <bitset>
#include <limits>
#include <algorithm>


using u64 = std::uint64_t;
using i64 = std::int64_t;
using f64 = double;


struct synthetic_trace {
  
  static constexpr u64 I = 4; // bytes per instruction
  static constexpr u64 block_avg = 7; // instructions per branch
  static constexpr f64 fwd_fraction = 0.9;
  static constexpr f64 loop_count_median = 8;
  static constexpr f64 loop_count_mode = 4;
  static constexpr f64 fwd_taken_proba = 0.3;
  static constexpr u64 region = 64; // region size (instructions) 
  static constexpr f64 diverge_fraction = 0.3;
  static_assert(block_avg >= 3);
  static_assert(loop_count_mode <= loop_count_median);

  const u64 footprint;
  const u64 nesting;
  
  static constexpr u64 NR = 3;
  std::mt19937_64 rng[NR];
  using bernoulli = std::bernoulli_distribution;
  bernoulli bias_dist {fwd_taken_proba};
  std::lognormal_distribution<f64> loop_count_dist {std::log(loop_count_median),std::sqrt(std::log(loop_count_median/loop_count_mode))};
  bernoulli noise_dist {0};

  u64 zone_address = 0;
  u64 bi = 0; // current block index
  u64 pc = 0;
  u64 stats_branch = 0;
  u64 stats_loop = 0;
  std::vector<u64> edge;

  struct block {
    u64 branch_addr = 0;
    u64 jump_target_addr = 0;
    u64 target_block = 0;
    bool is_loop = 0;
    u64 loop_count = 0;
    u64 loop_iter = 0;
    bool bias = 0;

    block(synthetic_trace &trace, u64 branchpc, u64 targetaddr)
    {
      branch_addr = branchpc;
      jump_target_addr = targetaddr;
      is_loop = jump_target_addr <= branch_addr;
      loop_count = std::ceil(trace.loop_count_dist(trace.rng[0]));
      bias = trace.bias_dist(trace.rng[1]);
    }

    bool loop()
    {
      assert(loop_count!=0);
      loop_iter++;
      if (loop_iter == loop_count) loop_iter = 0;
      return loop_iter != 0;
    }
  };

  std::vector<block> blocks;

  struct pattern {
    std::vector<bool> p;
    u64 hand = 0;

    pattern(u64 pattern_length)
    {
      p.resize(pattern_length);
    }

    void init(auto &rng)
    {
      bernoulli pattern_dist {diverge_fraction};
      for (auto &&e : p) e = pattern_dist(rng);
    }

    bool next()
    {
      hand++;
      if (hand==p.size()) hand=0;
      return p[hand];
    }
  };

  std::vector<pattern> patt;

  synthetic_trace(u64 footprint=1000, u64 patterns=8, u64 pattern_length=101, u64 nesting=5, u64 seed=0, f64 noise=0.001) : footprint(footprint), nesting(nesting), noise_dist{noise}, edge(footprint), patt(patterns,pattern_length)
  {
    assert(nesting * (block_avg-2) >= 2);
    u64 seeds[NR];
    std::seed_seq ss {seed};
    ss.generate(seeds,seeds+NR);
    for (u64 i=0; i<NR; i++) rng[i].seed(seeds[i]);
    zone_address = rng[0]();
    pc = zone_address;
    generate_edges(rng[0],rng[1]);
    for (u64 i=0; i<footprint-1; i++) {
      if (edge[i]==0) continue;
      blocks.emplace_back(*this,branch_address(i),branch_address(edge[i]));
    }
    // end the zone with a backward jump to the start
    blocks.emplace_back(*this,branch_address(footprint-1),zone_address);
    for (u64 i=0; i<blocks.size(); i++) {
      blocks[i].target_block = find_block(blocks[i].jump_target_addr);
    }
    assert(blocks.back().is_loop);
    blocks.back().loop_count = std::numeric_limits<u64>::max();
    for (auto &p : patt) p.init(rng[0]);
  }

  void generate_edges(auto &rng1, auto &rng2)
  {
    f64 step_proba = 2 * (1+1./nesting) / block_avg;
    assert(step_proba <= 1);
    bernoulli is_step {step_proba};
    bernoulli step_up;
    // set backward edges
    u64 level = 0;
    std::vector<u64> edge_end (nesting,0);
    for (u64 i=1; i<footprint; i++) {
      if (! is_step(rng1)) continue;
      u64 prev_level = level;
      if (step_up(rng2)) {
	if (level < nesting) level++;
      } else {
	if (level != 0) level--;
      }
      if (level < prev_level) edge[i] = edge_end.at(level); // set an edge
      if (level > prev_level) edge_end.at(prev_level) = i; // record an edge end
    }
    // revert some edges in proportion of forward jumps
    bernoulli is_fwd {fwd_fraction};
    for (u64 i=0; i<footprint; i++) {
      if (edge[i]!=0 && is_fwd(rng1)) {
	// revert the edge
	assert(edge[i]<i);
	assert(edge[edge[i]]==0);
	edge[edge[i]] = i;
	edge[i] = 0;
      }
    }
  }

  u64 branch_address(u64 offset) const
  {
    return zone_address + offset;
  }

  void print_stats(std::string s = "", std::ostream & os = std::cout) const
  {
    os << s << "branch sites: " << blocks.size() << std::endl;
    f64 loopsites = std::count_if(blocks.begin(),blocks.end(),[](const block &b){return b.is_loop;});
    os << std::setprecision(3);
    os << s << "loop sites fraction: " << loopsites/blocks.size() << std::endl;
    os << s << "loop dynamic fraction: " << f64(stats_loop)/stats_branch << std::endl;
  }

  u64 find_block(u64 address) const
  {
    // find the block containing the address
    assert(address >= zone_address);
    assert(address <= blocks.back().branch_addr);
    u64 low = 0;
    u64 high = blocks.size()-1;
    while ((low+1) < high) {
      u64 mid = (low+high)/2;
      if (address <= blocks[mid].branch_addr) {
	high = mid;
      } else {
	low = mid;
      }
    }
    return (address <= blocks[low].branch_addr)? low : high;
  }

  bool branch_dir(block &b)
  {
    if (b.is_loop) return b.loop();
    u64 patt_idx = (b.branch_addr/region) % std::size(patt);
    bool diverge = patt[patt_idx].next();
    bool random_flip = noise_dist(rng[0]);
    return b.bias ^ diverge ^ random_flip;
  }

  auto next_branch(bool byteaddr = true)
  {
    assert(bi < blocks.size());
    block& b = blocks[bi];
    stats_branch++;
    if (b.is_loop) stats_loop++;
    bool branch_taken = branch_dir(b);
    u64 next_addr = (branch_taken)? b.jump_target_addr : b.branch_addr + 1;
    bi = (branch_taken)? b.target_block : bi+1; // next block
    if (byteaddr) {
      return std::tuple{b.branch_addr*I, branch_taken, next_addr*I};
    } else {
      return std::tuple{b.branch_addr, branch_taken, next_addr};
    }
  }

  auto next_instruction()
  {
    u64 current_pc = pc;
    if (pc < blocks[bi].branch_addr) {
      pc++;
      return std::tuple{current_pc*I, /*branch*/false, /*taken*/false, /*next pc*/pc*I};
    } else {
      // branch
      auto [branch_pc/*unused*/,branch_taken,next_pc] = next_branch(false);
      pc = next_pc;
      return std::tuple{current_pc*I, /*branch*/true, branch_taken, /*next pc*/pc*I};
    }
  }
};

#endif // SYNTHETIC_TRACE_H
