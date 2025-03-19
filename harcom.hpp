// HARCOM: Hardware Complexity Model
// Language: C++20

// Copyright (c) 2025 INRIA
// Copyright (c) 2025 Pierre Michaud

// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge, publish, distribute,
// sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE
// AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#ifndef HARCOM_H
#define HARCOM_H

#include <cstdint>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <ios>
#include <type_traits>
#include <concepts>
#include <bitset>
#include <bit>
#include <string>
#include <utility>
#include <initializer_list>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <array>
#include <vector>
#include <limits>


class simulator; // has access to private members of class "val"


namespace hcm {

  using u8 = std::uint8_t;
  using u16 = std::uint16_t;
  using u32 = std::uint32_t;
  using u64 = std::uint64_t;
  using i8 = std::int8_t;
  using i16 = std::int16_t;
  using i32 = std::int32_t;
  using i64 = std::int64_t;
  using f32 = float;
  using f64 = double;
  
  template<typename T>
  concept arith = std::integral<T> || std::floating_point<T>;

  template<typename T>
  concept action = requires (T f) {f();} || requires (T f, u64 i) {f(i);};

  template<u64 N, arith T> class val;
  template<u64 N, arith T> class reg;
  
  template<typename T>
  concept valtype = requires(const T &x) {[]<u64 N, arith U>(const val<N,U>&){}(x);};

  template<valtype T, u64 N> class arr;

  template<typename T>
  concept arrtype = requires(const T &x) {[]<valtype U, u64 N>(const arr<U,N>&){}(x);};

  template<typename T>
  concept regtype = requires(const T &x) {[]<u64 N, arith U>(const reg<N,U>&){}(x);};
  
  template<typename T>
  concept memdatatype = (valtype<T> || arrtype<T>) && ! regtype<T> && ! regtype<typename T::type>;

  
  constexpr auto min(auto a, auto b)
  {
    return (b<a)? b : a;
  }

  constexpr auto max(auto a, auto b)
  {
    return (b>a)? b : a;
  }  

  auto to_unsigned(std::integral auto x) {return std::make_unsigned_t<decltype(x)>(x);}
  auto to_unsigned(f32 x) {return std::bit_cast<u32>(x);}
  auto to_unsigned(f64 x) {return std::bit_cast<u64>(x);}

  template<u64 N, std::integral T>
  auto truncate(T x)
  {
    if constexpr (N >= sizeof(T)*8) {
      return x;
    } else {
      return x & ((std::make_unsigned_t<T>(1)<<N)-1);
    }
  }

  template<u64 N, std::integral T>
  auto ones(T x)
  {
    return std::popcount(truncate<N>(to_unsigned(x)));
  }
  
  constexpr u64 to_pow2(f64 x)
  {
    // returns the power of 2 closest to x
    u64 n0 = std::bit_floor(u64(llround(floor(x))));
    u64 n1 = std::bit_ceil(u64(llround(ceil(x))));
    return (x*x > n0*n1)? n1 : n0;
  }

  constexpr u8 bit_reversal(u8 x)
  {
    u8 y = x&1;
    for (int i=1; i<8; i++) {
      x >>= 1;
      y = (y<<1) | x&1;
    }
    return y;
  }

  static constexpr auto reversed_byte = [] () {
    std::array<u8,256> b;
    for (u64 i=0; i<256; i++) {
      b[i] = bit_reversal(i);
    }
    return b;
  }();

  auto constexpr reverse_bits(std::unsigned_integral auto x)
  {
    using T = decltype(x);
    T y = reversed_byte[x & 0xFF];
    for (int i=1; i<sizeof(T); i++) {
      x >>= 8;
      y = (y<<8) | reversed_byte[x & 0xFF];
    }
    return y;
  }


  // ###########################

  // utility for doing a static loop over an integer template argument (0,...,N-1)
  // the loop body is a C++ lambda with one integer template parameter
  
  template<typename T>
  concept static_loop_body = requires (T obj) {obj.template operator()<0>();};

  template <u64 ...I>
  constexpr void static_loop(static_loop_body auto F, std::integer_sequence<u64,I...>)
  {
    (F.template operator()<I>(),...);
  }

  template<u64 N>
  constexpr void static_loop(static_loop_body auto F)
  {
    static_loop(F,std::make_integer_sequence<u64,N>{});
  }
  

  // ###########################

  // wires: see Nikolic, FPGA 2021
  // assume unique linear capacitance (fF/um) for all metal layers
  inline constexpr f64 METALCAP_fF = 0.2; // fF/um
  inline constexpr f64 METALRES[] = {150/*Mx*/,25/*My*/}; // Ohm/um
  
  inline constexpr f64 PINV = 1; // ratio of drain capacitance to gate capacitance

  // TSMC 5nm, Chang, IEEE JSSC, jan 2021: 
  // each transistor in 6T cell has single fin (for high density)
  // bitline capacitance per cell = 1 drain capacitance (Cd) + 1 Mx wire segment (Cw)
  // Flying bitline (Mx+2) reduces bitline capacitance: 50%(bottom)/65%(top)
  // bottom bitline = (N/2)x(Cd+Cw)
  // top bitline = (N/2)x(Cd+Cw) + (N/2)xCw
  // (Cd+2Cw)/(Cd+Cw) = 65/50 ==> Cd/Cw = 0.7/0.3 = 2.33
  inline constexpr struct {
    u64 transistors = 6; // single RW port
    f64 area = 0.02; // um^2
    f64 density = transistors / area;
    f64 aspect_ratio = 2; // wordline longer than bitline
    f64 drain_blwire_cap_ratio = 2.33;
    f64 bitline_length = sqrt(area/aspect_ratio); // um
    f64 wordline_length = bitline_length * aspect_ratio; // um
    f64 bitline_wirecap = bitline_length * METALCAP_fF; // fF
    f64 wordline_wirecap = wordline_length * METALCAP_fF; // fF
    f64 bitline_resistance = bitline_length * METALRES[0]; // Ohm
    f64 wordline_resistance = wordline_length * METALRES[0]; // Ohm
    f64 drain_capacitance = drain_blwire_cap_ratio * bitline_wirecap; // fF
    f64 gate_capacitance = drain_capacitance / PINV; // fF
    f64 bitline_capacitance = bitline_wirecap + drain_capacitance; // fF
    f64 wordline_capacitance = wordline_wirecap + 2 * gate_capacitance; // fF

    void print(std::ostream & os = std::cout) const
    {
      os << "gate capacitance (fF): " << gate_capacitance << std::endl;
      os << "drain capacitance (fF): " << drain_capacitance << std::endl;
      os << "wordline (per cell):" << std::endl;
      os << "  length (um): " << wordline_length << std::endl;
      os << "  wire resistance (Ohm): " << wordline_resistance << std::endl;
      os << "  wire capacitance (fF): " << wordline_wirecap << std::endl;
      os << "  capacitance (fF): " << wordline_capacitance << std::endl;
      os << "bitline (per cell):" << std::endl;
      os << "  length (um): " << bitline_length << std::endl;
      os << "  wire resistance (Ohm): " << bitline_resistance << std::endl;
      os << "  wire capacitance (fF): " << bitline_wirecap << std::endl;
      os << "  capacitance (fF): " << bitline_capacitance << std::endl;      
    }
  } SRAM_CELL;

  // CGATE = gate capacitance of single-fin NMOS
  inline constexpr f64 CGATE_fF = SRAM_CELL.gate_capacitance; // fF
  inline constexpr f64 CGATE_pF = CGATE_fF * 1e-3; // pF
  inline constexpr f64 METALCAP = METALCAP_fF / CGATE_fF; // linear metal capacitance in units of CGATE
  inline constexpr f64 METALCAP_pF = METALCAP * CGATE_pF; // pF
  inline constexpr f64 VDD = 0.75; // supply voltage (V)
  inline constexpr f64 IDSAT = 60e-6; // NMOS Idsat per fin (A) (IRDS 2021, Ioff ~= 10 pA/fin)
  // IEFF, REFF ==> Razavieh, Device Research Conference, 2018
  inline constexpr f64 IEFF = IDSAT/2; // NMOS effective current per fin (A)
  inline constexpr f64 REFF = VDD/(2*IEFF); // NMOS effective resistance (Ohm)
  inline constexpr f64 GAMMA = 1; // ratio of PMOS to NMOS capacitance at same conductance
  inline constexpr f64 INVCAP = 1+GAMMA; // input capacitance of single-fin inverter relative to CGATE
  inline constexpr f64 TAU_ps = CGATE_pF * REFF; // intrinsic delay (ps)

  inline constexpr f64 RANDOM_LOGIC_DENSITY = 0.25 * SRAM_CELL.density;
  
  // Currently, we assume single-fin transistors everywhere (except for large fanout).
  // TODO: it would be more realistic to assume 2-fin transistors (except in SRAM cell).
  // As we do not model the effect of interconnects on gates delay, this does not matter for delay,
  // however we underestimate the energy (2-fin transistor has larger gate capacitance)

  inline constexpr f64 DSE = 6; // default stage effort (delay vs energy tradeoff)
  inline constexpr u64 DSMAX = 100; // default maximum scale (FIXME)
  
  
  constexpr f64 energy_fJ(f64 cap_fF, f64 vdiff/*volt*/)
  {
    // energy dissipated for changing capacitance voltage by VDIFF
    assert(vdiff>0);
    return cap_fF * vdiff * (VDD - vdiff/2); // fJ
  }


  constexpr f64 wire_res_delay(f64 res/*Ohm*/, f64 wirecap_pF, f64 loadcap_pF=0)
  {
    // approximate distributed RC line as lumped PI1 circuit (Rao, DAC 1995)
    f64 ElmoreDelay = res * (wirecap_pF/2 + loadcap_pF); // ps
    // RC delay approximately equal to Elmore delay x ln(2) (Gupta et al., IEEE TCAD jan 1997)
    return ElmoreDelay * log(2); // ps
  }


  constexpr f64 proba_switch(f64 bias)
  {
    assert(bias>=0 && bias<=1);
    return 2 * bias * (1-bias);
  }


  constexpr f64 proba_bias(f64 proba_switch)
  {
    assert(proba_switch <= 0.5); // random switching (Bernoulli process)
    return 0.5 * (1-sqrt(1-2*proba_switch)); // does not exceed 1/2
  }
  
  
  struct circuit {
    u64 t = 0; // transistors
    f64 d = 0; // delay (ps)
    f64 ci = 0; // maximum input capacitance relative to CGATE
    f64 cg = 0; // total gate capacitance (all transistors) relative to CGATE
    f64 e = 0; // energy (fJ)
    f64 w = 0; // total wiring (um)

    constexpr circuit() : t(0), d(0), ci(0), cg(0), e(0), w(0) {}
    
    constexpr circuit(u64 t, f64 d, f64 ci, f64 cg, f64 bias) : t(t), d(d), ci(ci), cg(cg), e(gate_energy(bias)), w(0)
    {
      // for gates only
      assert(t!=0);
      assert(d>0);
      assert(ci>0);
      assert(cg>0);
      assert(e>=0);
    }

    constexpr bool nogate() const
    {
      return t==0;
    }

    u64 delay() const {return llround(d);}
    f64 energy() const {return e;}
    
    void print(std::string s = "", std::ostream & os = std::cout) const
    {
      os << s;
      os << std::setprecision(4);
      os << "xtors: " << t;
      os << " ; wires (um): " << w;
      os << " ; ps: " << d;
      os << " ; fJ: " << e;
      os << " ; input cap: " << ci;
      //os << " ; gate cap: " << cg;
      os << std::endl;
    }

    constexpr f64 gate_energy(f64 bias) const
    {
      // TODO: static power, short circuit
      f64 c = (1+PINV) * CGATE_fF * cg; // switching capacitance (fF)
      // FIXME: not all drain capacitances switch
      return energy_fJ(c,VDD) * proba_switch(bias); // fJ
    }

    constexpr f64 cost() const
    {
      return e*d*d;
    }    
    
    constexpr circuit operator+ (const circuit &x) const
    {
      // series (LHS output connected to RHS input)
      circuit chain;
      chain.t = t + x.t;
      chain.d = d + x.d;
      chain.ci = (ci==0)? x.ci : ci;
      chain.cg = cg + x.cg;
      chain.e = e + x.e;
      chain.w = w + x.w;
      return chain;
    }

    constexpr circuit operator* (f64 n) const
    {
      // parallel, distinct inputs
      if (n==0) return {};
      circuit rep;
      rep.t = t * n;
      rep.d = d;
      rep.ci = ci;
      rep.cg = cg * n;
      rep.e = e * n;
      rep.w = w * n;
      return rep;
    }
    
    constexpr circuit operator|| (const circuit &x) const
    {
      // parallel, distinct inputs
      circuit para;
      para.t = t + x.t;
      para.d = max(d,x.d);
      para.ci = max(ci,x.ci);
      para.cg = cg + x.cg;
      para.e = e + x.e;
      para.w = w + x.w;
      return para;
    }
    
    constexpr circuit operator| (const circuit &x) const
    {
      // parallel, single input
      circuit para;
      para.t = t + x.t;
      para.d = max(d,x.d);
      para.ci = ci + x.ci;
      para.cg = cg + x.cg;
      para.e = e + x.e;
      para.w = w + x.w;
      return para;
    }
  };

  
  struct basic_gate {
    using cilist = std::initializer_list<f64>;
    static constexpr u64 MAX_INPUT_TYPES = 3;
    u64 tr = 0; // number of transistors
    f64 ci[MAX_INPUT_TYPES] = {0}; // inputs capacitances, relative to CGATE
    f64 cp = 0; // parasitic output capacitance, relative to CGATE
    f64 cg = 0; // total gate capacitances, relative to CGATE
    u64 nci = 0; // number of input types
    f64 cimax = 0; // maximum input capacitance
    
    constexpr void init_ci(cilist l)
    {
      nci = l.size();
      assert(nci!=0);
      assert(nci<=MAX_INPUT_TYPES);
      cimax = std::max(l);
      int i = 0;
      for (auto e : l) ci[i++] = e;
    }

    constexpr basic_gate() {}
    
    constexpr basic_gate(u64 tr, cilist cl, f64 p, f64 cg) : tr(tr), cp(p*PINV), cg(cg)
    {
      init_ci(cl);
    }

    template<u64 I=0>
    constexpr f64 icap(f64 scale=1) const
    {
      assert(I<nci);
      return ci[I] * scale;
    }

    template<u64 I=0>
    constexpr f64 logical_effort() const
    {
      assert(I<nci);
      return ci[I] / INVCAP;
    }
    
    constexpr circuit make(f64 co, f64 scale=1, f64 bias=0.5) const
    {
      //  co = output (load) capacitance relative to CGATE
      assert(scale>=1);
      assert(bias>=0 && bias<=1);
      assert(co>=0);
      f64 delay_ps = (cp + co/scale) * TAU_ps;
      return {tr, delay_ps, cimax*scale, cg*scale, bias};
    }
  };

  
  struct inv : basic_gate { // inverter
    constexpr inv() : basic_gate(2,{INVCAP},INVCAP,INVCAP) {}
  };

  struct nand : basic_gate { // single-input NAND is inverter
    constexpr nand(u64 n) : basic_gate(2*n,{n+GAMMA},n+n*GAMMA,n*(n+GAMMA)) {assert(n>=1);}
  };

  struct nor : basic_gate { // single-input NOR is inverter
    constexpr nor(u64 n) : basic_gate(2*n,{1+n*GAMMA},n+n*GAMMA,n*(1+n*GAMMA)) {assert(n>=1);}
  };

  struct and_nor /*aka AOI21*/: basic_gate { // ~(a|(b&c))
    constexpr and_nor() : basic_gate(6,{1+2*GAMMA/*a*/,2+2*GAMMA/*b*/,2+2*GAMMA/*c*/},3+2*GAMMA,5+6*GAMMA) {}
  };

  struct or_nand /*aka OAI21*/: basic_gate { // ~(a&(b|c))
    constexpr or_nand() : basic_gate(6,{2+GAMMA/*a*/,2+2*GAMMA/*b*/,2+2*GAMMA/*c*/},2+3*GAMMA,6+5*GAMMA) {}
  };

  struct mux_inv_tri : basic_gate { // inverting MUX (tristate inverters)
    constexpr mux_inv_tri(u64 n) : basic_gate(4*n,{2+2*GAMMA/*data*/,2/*sel*/,2*GAMMA/*csel*/},n*(2+2*GAMMA),n*(4+4*GAMMA)) {assert(n>=1);}
  };

  struct inv_tri : mux_inv_tri { // tristate inverter
    constexpr inv_tri() : mux_inv_tri(1) {}
  };
  
  struct xor_cpl : basic_gate { // a,~a,b,~b ==> a^b
    constexpr xor_cpl() : basic_gate(8,{2+2*GAMMA},4+4*GAMMA,8+8*GAMMA) {}
  };

  using xnor_cpl = xor_cpl; // a,~a,b,~b ==> ~(a^b)
  

  template<f64 STAGE_EFFORT=DSE>
  constexpr u64 num_stages(f64 path_effort, bool odd=false)
  {
    // higher stage effort increases delay but reduces energy
    assert(path_effort>=1);
    u64 n = log(path_effort) / log(STAGE_EFFORT);
    assert(n>=0);
    n += (n & 1) ^ odd;
    return n;
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit buffer(f64 co, bool cpl, f64 scale=1, f64 bias=0.5)
  {
    f64 ci = inv{}.icap(scale);
    f64 fo = min(max(1,co/ci),SMAX);
    u64 ninv = num_stages<SE>(fo,cpl);
    circuit buf;
    if (ninv!=0) {
      f64 foi = pow(fo,1./ninv);
      assert(foi>=1);
      for (int i=0; i<ninv; i++) {
	ci *= foi;
	buf = buf + inv{}.make((i==(ninv-1))? co:ci,scale,bias);
	scale *= foi;
      }
    }
    if (buf.nogate()) buf.ci = co;
    return buf;
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX, u64 L=1/*layer*/>
  constexpr circuit wire(f64 length/*um*/, bool cpl=0, f64 cload=0, bool dload=true, f64 bias=0.5)
  {
    // Unidirectional segmented wire
    // The wire is divided into equal length segments
    // The first segment is driven by a buffer, other segments by a large inverter
    // Load capacitance can be terminal (dload=0) or uniformly distributed (dload=1)
    // TODO: bidirectional wire is possible (with tristate inverters), but can it be pipelined?
    static_assert(L<std::size(METALRES));
    assert(length>=0);
    if (length==0) dload = false;
    f64 co = (dload)? 0 : cload; // terminal load capacitance, relative to CGATE
    // for distributed load cap, merge with wire capacitance (FIXME: is this good approx?)
    f64 linearcap = METALCAP + ((dload)? (cload/length) : 0); // CGATE/um
    f64 seglen = sqrt(2*REFF*INVCAP*(1+PINV) / (linearcap*METALRES[L])); // um
    f64 optinvscale = sqrt(REFF * linearcap / (INVCAP * METALRES[L]));
    f64 invscale = min(SMAX,optinvscale); // inverter scale
    u64 nseg = max(1, llround(length / seglen)); // number of segments
    seglen = length / nseg; // segment length (um)
    f64 cseg = linearcap * seglen; // wire segment capacitance relative to CGATE
    // transistors:
    f64 cfirst = cseg + ((nseg==1)? co : 0);
    circuit c = buffer<SE,SMAX>(cfirst,(nseg-1+cpl)&1,1,bias);
    circuit rep = inv{}.make(cseg,invscale,bias); // repeater = inverter
    for (int i=1; i<(nseg-1); i++) c = c + rep;
    if (nseg>1) c = c + inv{}.make(cseg+co,invscale,bias); // last segment drives output cap
    // wire:
    c.w = length;
    c.d += wire_res_delay(METALRES[L]*seglen, cseg*CGATE_pF, rep.ci*CGATE_pF) * (nseg-1);
    c.d += wire_res_delay(METALRES[L]*seglen, cseg*CGATE_pF, co*CGATE_pF); // last segment
    c.e += energy_fJ(METALCAP_fF*length,VDD) * proba_switch(bias);
    if (c.nogate()) c.ci = cfirst; // just a wire
    return c;
  }

  
  constexpr circuit majority_nand3(f64 co, f64 scale=1, f64 bias=0.5)
  {
    // a,b,c ==> ab+ac+bc = ~(~(ab)~(ac)~(bc))
    circuit na3 = nand{3}.make(co,scale,bias);
    circuit na = nand{2}.make(na3.ci,scale,bias);
    return ((na | na) || na) + na3;
  }

  
  constexpr circuit majority_aoi21(f64 co, f64 scale=1, f64 bias=0.5)
  {
    // a,b,c ==> ab+ac+bc = ~(~(b+c)|(~a&~(bc)))
    and_nor gate;
    f64 c1 = gate.icap<0>(scale);
    f64 c2 = gate.icap<1>(scale);
    circuit i = inv{}.make(c2,scale,bias);
    circuit na = nand{2}.make(c2,scale,bias);
    circuit no = nor{2}.make(c1,scale,bias);
    circuit maj = ((na | no) || i) + gate.make(co,scale,bias); // TODO: and_nor inputs bias
    return maj;
  }

  
  constexpr circuit majority(f64 co, f64 scale=1, f64 bias=0.5)
  {
    return majority_aoi21(co,scale,bias);
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX, u64 ARITY=4>
  constexpr circuit nand_nor_tree(u64 n, bool nandfirst, bool cpl, f64 co, f64 scale=1, f64 bias=0.5)
  {
    // alternate NANDs and NORs
    static_assert(ARITY >= 2);
    if (n==0) return {};
    if (n==1) return (cpl)? inv{}.make(co,scale,bias) : circuit{};
    basic_gate gate[2][ARITY+1];
    for (int i=1; i<=ARITY; i++) {
      gate[0][i] = nor(i);
      gate[1][i] = nand(i);
    }
    u64 ngates[2][ARITY+1];

    auto partition = [] (u64 n) {
      u64 w = min(n,ARITY);
      u64 r = n%w;
      u64 nw = n/w;
      if (r!=0) {
	r += w;
	nw--;
      }
      return std::tuple{w,nw,r-r/2,r/2};
    };
    
    auto populate_gates = [partition,&ngates] (bool nand_stage, u64 n) {
      assert(n>0);
      if (n==1) return u64(0);
      auto [w,nw,w2,w3] = partition(n);
      assert(w<=ARITY && w2<=ARITY && w3<=ARITY);
      for (u64 &e : ngates[nand_stage]) e = 0;
      ngates[nand_stage][w] += nw;
      ngates[nand_stage][w2]++;
      ngates[nand_stage][w3]++;
      return (nw==0)? w2 : w; // widest
    };

    auto path_logical_effort = [nandfirst,gate,partition] (u64 n) {
      f64 path_le = 1;
      bool nand_stage = nandfirst;
      while (n>1) {
	auto [w,nw,w2,w3] = partition(n);
	u64 width = (nw==0)? w2 : w;
	assert(width!=0);
	path_le *= gate[nand_stage][width].logical_effort();
	n = nw;
	if (w2) n++;
	if (w3) n++;
	nand_stage ^= 1;
      }
      return path_le;
    };

    bool nand_stage = nandfirst;
    u64 width = populate_gates(nand_stage,n);
    f64 ci = gate[nand_stage][width].icap(scale);
    f64 fanout = min(co/ci,SMAX);
    u64 path_effort = max(1,path_logical_effort(n)*fanout);
    u64 depth_target = num_stages<SE>(path_effort);
    u64 depth = llround(ceil(log(n)/log(ARITY)));
    bool extra_inv = (depth & 1) ^ cpl; // odd number of stages if cpl=true, even otherwise
    depth += extra_inv;
    u64 ninv = 0;
    if (depth < depth_target) {
      // add inverters
      ninv = depth_target - depth;
      ninv += ninv & 1; // even
      depth += ninv;
    }
    ninv += extra_inv;
    assert(depth!=0);
    f64 stage_effort = pow(path_effort,1./depth);

    // build the tree
    circuit tree;
    n = std::accumulate(ngates[nand_stage]+1,ngates[nand_stage]+1+ARITY,0);
    do {
      assert(width>=2);
      u64 prev_width = width;
      width = populate_gates(nand_stage^1,n);
      f64 next_le = (width!=0)? gate[nand_stage^1][width].logical_effort() : 1;
      f64 next_scale = max(1,scale * stage_effort / next_le);
      f64 cload = (width!=0)? gate[nand_stage^1][width].icap(next_scale) : (depth>1)? inv{}.icap(next_scale) : co;
      circuit stage;
      for (u64 w=1; w<=ARITY; w++) {
	stage = stage || (gate[nand_stage][w].make(cload,scale,bias) * ngates[nand_stage][w]);
      }
      tree = tree + stage;
      bias = pow(bias,prev_width); // each NAND/NOR stage reduces switching probability
      depth--;
      scale = next_scale;
      nand_stage ^= 1;
      if (width!=0) n = std::accumulate(ngates[nand_stage]+1,ngates[nand_stage]+1+ARITY,0);
    } while (width!=0);
    
    if (ninv!=0) {
      // add inverters
      for (int i=0; i<ninv; i++) {
	f64 next_scale = scale * stage_effort;
	f64 cload = (depth>1)? inv{}.icap(next_scale) : co;
	tree = tree + inv{}.make(cload,scale,bias);
	depth--;
	scale = next_scale;
      }
    }
    assert(depth==0);
    return tree;
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit anding(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    return nand_nor_tree<SE,SMAX>(n,true,false,co,scale,bias);
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit oring(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    return nand_nor_tree<SE,SMAX>(n,false,false,co,scale,bias);
  }  

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit nanding(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    return nand_nor_tree<SE,SMAX>(n,true,true,co,scale,bias);
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit noring(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    return nand_nor_tree<SE,SMAX>(n,false,true,co,scale,bias);
  }  

  
  constexpr circuit xor2(f64 co, f64 scale=1, f64 bias=0.5)
  {
    circuit x = xor_cpl{}.make(co,scale,bias);
    circuit i = inv{}.make(x.ci,scale,bias);
    circuit c = i*2+x;
    c.ci = i.ci + x.ci;
    return c;
  }

  
  constexpr circuit xnor2(f64 co, f64 scale=1, f64 bias=0.5)
  {
    return xor2(co,scale,bias);
  }  

  
  constexpr circuit parity_small(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    // XOR tree
    if (n<=1) return {};
    f64 ci = xor2(1/*whatever*/,scale).ci;
    circuit tree;
    while (n>1) {
      if (n==2) {
	// last stage
	tree = tree + xor2(co,scale,bias);
      } else {
	// not the last stage
	tree = tree + xor2(ci,scale,bias) * (n/2);
      }
      n -= n/2;
      bias = 2*bias*(1-bias);
    }
    return tree;
  }

  
  constexpr circuit parity_fast(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    // XOR/XNOR tree
    if (n<=1) return {};
    f64 ci = xor_cpl{}.icap(scale) + xnor_cpl{}.icap(scale);
    circuit tree = inv{}.make(ci,scale,bias) * n;
    tree.ci += ci;
    while (n>1) {
      if (n==2) {
	// last stage
	tree = tree + xor_cpl{}.make(co,scale,bias);
      } else {
	tree = tree + (xor_cpl{}.make(ci,scale,bias) | xnor_cpl{}.make(ci,scale,bias)) * (n/2);
      }
      n -= n/2;
      bias = 2*bias*(1-bias);
    }
    return tree;
  }

  
  constexpr circuit parity(u64 n, f64 co, f64 scale=1, f64 bias=0.5)
  {
    return parity_small(n,co,scale,bias);
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit decode1(u64 no/*outputs*/, f64 co, f64 pitch=0/*um*/)
  {
    // single-level decoder
    // co = load capacitance per output (relative to CGATE)
    // pitch = distance between outputs (um)
    assert(no!=0);
    if (no==1) return buffer<SE,SMAX>(co,false);
    assert(no>=2);
    u64 ni = std::bit_width(no-1); // inputs (address bits)
    if (ni == 1) {
      return wire<SE,SMAX>(pitch/2,false,co,false) | wire<SE,SMAX>(pitch/2,true,co,false);
    }
    circuit g = anding<SE,SMAX>(ni,co);
    f64 h = pitch * (no-1); // decoder height (um)
    // each input is connected to no/2 gates
    circuit wireup = wire<SE,SMAX>(h/2,false,g.ci*(no/4)) | wire<SE,SMAX>(h/2,true,g.ci*(no/4));
    circuit wiredown = wire<SE,SMAX>(h/2,false,g.ci*(no/2-no/4)) | wire<SE,SMAX>(h/2,true,g.ci*(no/2-no/4));
    circuit dec = (wireup | wiredown) * ni + g * no;
    return dec;
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit decode2(u64 no/*outputs*/, f64 co, f64 pitch=0/*um*/)
  {
    // two-level decoder
    // co = load capacitance per output (relative to CGATE)
    // pitch = distance between outputs (um)
    if (no<4) return decode1<SE,SMAX>(no,co,pitch);
    u64 ni = std::bit_width(no-1); // inputs (address bits)
    // 2 predecoders
    u64 pni[2] = {ni/2, ni-ni/2}; // address bits per predecoder
    u64 pno[2] = {1ull<<pni[0], 1ull<<pni[1]}; // ouputs per predecoder
    f64 bias[2] = {1./pno[0], 1./pno[1]}; // one output activated per predecoder
    circuit g = anding<SE,SMAX>(2,co,1,(bias[0]+bias[1])/2);
    f64 h = pitch * (no-1); // decoder height (um)
    circuit pdwire[2] = {wire<SE,SMAX>(h/2,false,g.ci*(pno[1]/2),true,bias[0]) | wire<SE,SMAX>(h/2,false,g.ci*(pno[1]-pno[1]/2),true,bias[0]), wire<SE,SMAX>(h/2,false,g.ci*(pno[0]/2),true,bias[1]) | wire<SE,SMAX>(h/2,false,g.ci*(pno[0]-pno[0]/2),true,bias[1])};
    circuit pdec[2] = {decode1<SE,SMAX>(pno[0],pdwire[0].ci), decode1<SE,SMAX>(pno[1],pdwire[1].ci)};
    return ((pdec[0]+(pdwire[0]*pno[0])) || (pdec[1]+(pdwire[1]*pno[1]))) + g * no;
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit decode2_rep(u64 no/*outputs*/, f64 co, f64 pitch=0/*um*/, u64 rep=1)
  {
    // two-level decoder with replicated outputs
    // rep = number of times an output is replicated
    // pitch = distance between outputs/replicas (um)
    assert(no>=4);
    u64 ni = std::bit_width(no-1); // inputs (address bits)
    // 2 predecoders
    u64 pni[2] = {ni/2, ni-ni/2}; // address bits per predecoder
    u64 pno[2] = {1ull<<pni[0], 1ull<<pni[1]}; // ouputs per predecoder
    f64 bias[2] = {1./pno[0], 1./pno[1]}; // one output activated per predecoder
    circuit g = anding<SE,SMAX>(2,co,1,(bias[0]+bias[1])/2);
    f64 h = pitch * (no*rep-1); // decoder height (um)
    circuit pdwire[2] = {wire<SE,SMAX>(h/2,false,g.ci*(pno[1]/2)*rep,true,bias[0]) | wire<SE,SMAX>(h/2,false,g.ci*(pno[1]-pno[1]/2)*rep,true,bias[0]), wire<SE,SMAX>(h/2,false,g.ci*(pno[0]/2)*rep,true,bias[1]) | wire<SE,SMAX>(h/2,false,g.ci*(pno[0]-pno[0]/2)*rep,true,bias[1])};
    circuit pdec[2] = {decode1<SE,SMAX>(pno[0],pdwire[0].ci), decode1<SE,SMAX>(pno[1],pdwire[1].ci)};
    return ((pdec[0]+(pdwire[0]*pno[0])) || (pdec[1]+(pdwire[1]*pno[1]))) + g * (no*rep);
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX, u64 ARITY=4>
  constexpr auto mux_tree(u64 n/*inputs (data)*/, u64 m/*bits per input*/, f64 co)
  {
    // select input in encoded form
    // made from inverting MUXes ==> number of stages must be even
    // FIXME: distances and wires are not modeled here
    static_assert(ARITY>=2 && std::has_single_bit(ARITY));
    assert(n>=2);
    assert(m>=1);

    auto multiplexer = [] (u64 w, f64 cload) {
      assert(w>0);
      return (w==1)? inv{}.make(cload) : mux_inv_tri{w}.make(cload);
    };
    
    auto selector = [m] (u64 w, u64 nmuxes) {
      assert(w>1);
      mux_inv_tri muxi{w};
      if (w==2) {
	// single select signal
	f64 selcap = (muxi.icap<1>() + muxi.icap<2>()) * nmuxes * m;
	return buffer<SE,SMAX>(selcap,true) |  buffer<SE,SMAX>(selcap,false);
      } else {
	assert(w>2);
	f64 selcap = muxi.icap<1>() * nmuxes * m;
	f64 cselcap = muxi.icap<2>() * nmuxes * m;
	circuit sel = buffer<SE,SMAX>(selcap,true) | buffer<SE,SMAX>(selcap,false);
	return decode1<SE,SMAX>(w,sel.ci) + sel * w;
      }
    };

    f64 muxdatacap = mux_inv_tri{2}.icap<0>(); // independent of MUX width
    circuit sel;
    circuit data;
    u64 nstages = 0;

    while (n>1) {
      nstages++;
      if constexpr (ARITY>2) {
	if ((nstages & 1)==0 && (n+1)==ARITY) {
	  n++; // last stage = MUX with unused input
	}
      }
      u64 w = ARITY;
      while (w>n) w>>=1;
      u64 nw = n/w;
      u64 r = n%w;
      assert(nw>=1);
      bool last_mux = (nw==1) && (r==0);
      bool last_stage = last_mux && (nstages & 1)==0;
      f64 cload = (last_stage)? co : (last_mux)? inv{}.icap() : muxdatacap;
      circuit sel_stage;
      circuit mux_stage;
      if (std::has_single_bit(r)) {
	// one smaller MUX
	if (r==1) {
	  // just an inverter
	  mux_stage = mux_stage || inv{}.make(cload);
	} else {
	  sel_stage = sel_stage | selector(r,1);
	  mux_stage = mux_stage || multiplexer(r,cload);
	}
      } else if (r!=0) {
	// one extra MUX with unused inputs
	nw++;
	r = 0;
      }
      sel_stage = sel_stage | selector(w,nw);
      mux_stage = mux_stage || multiplexer(w,cload) * nw;
      mux_stage = mux_stage * m;
      sel = sel || sel_stage;
      data = data + mux_stage;
      n = nw + ((r!=0)? 1:0);
    }
    if (nstages & 1) {
      data = data + inv{}.make(co);
      nstages++;
    }
    return std::array{sel,data};
  }

  
  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr auto mux(u64 n/*inputs (data)*/, u64 m/*bits per input*/, f64 co)
  {
    assert(n>=2);
    assert(m>=1);
    auto mux2 = mux_tree<SE,SMAX,2>(n,m,co);
    auto mux4 = mux_tree<SE,SMAX,4>(n,m,co);
    return ((mux4[0]+mux4[1]).cost() < (mux2[0]+mux2[1]).cost())? mux4 : mux2;
  }


  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit grid_demux(u64 nbits, u64 nx, u64 ny, f64 dx/*um*/, f64 dy/*um*/, f64 cnode=0, u64 abits=1)
  {
    // sends an nbits-bit payload to one of nx*ny nodes
    // the nx*ny nodes form a regular grid
    // distance between nodes is dx (horizontal) and dy (vertical)
    // cnode = capacitance (per bit) of the selected node, relative to CGATE
    // binary tree network; 2-way forks implemented with tristate inverters for reducing wire energy
    // the packet includes the payload and the address bits to route itself through the network
    // when the packet passes through a fork, it drops one address bit (the fork-select one)
    // wires that are off the route should not switch
    // TODO: can we prevent glitching with artificial delays, or do we need clocked logic?
    assert(nbits!=0);
    assert(nx!=0 && ny!=0);
    assert(dx>=0 && dy>=0);
    assert(abits!=0);
    if (nx==1 && ny==1) {
      return {};
    }
    u64 packetbits = nbits + abits-1;
    f64 psel = 1./(nx*ny); // a single branch is selected
    assert(psel<=0.5);
    f64 bias = proba_bias(psel); // pseudo bias
    inv_tri itri; // tristate inverter
    f64 selcap = (itri.icap<1>() + itri.icap<2>()) * packetbits; // FIXME: wire not modeled
    circuit select = buffer<SE,SMAX>(selcap,false,1,bias) | buffer<SE,SMAX>(selcap,true,1,bias);
    if ((nx>=2 && dx<=dy) || ny==1) {
      // demux horizontally
      auto packwire = [=](bool cpl) {return wire(dx/2,cpl,cnode,false,proba_bias(psel*0.5));};
      circuit branch = (itri.make(packwire(true).ci,1,bias) + packwire(true)) * packetbits;
      circuit demux = (select + (branch | branch)) * ((nx/2)*ny);
      demux = demux || packwire(false) * ((nx&1) * packetbits);
      return grid_demux(nbits,(nx+1)/2,ny,dx*2,dy,demux.ci,abits+1) + demux;
    } else {
      // demux vertically
      assert(ny>=2);
      auto packwire = [=](bool cpl) {return wire(dy/2,cpl,cnode,false,proba_bias(psel*0.5));};
      circuit branch = (itri.make(packwire(true).ci,1,bias) + packwire(true)) * packetbits;
      circuit demux = (select + (branch | branch)) * (nx*(ny/2));
      demux = demux || packwire(false) * ((ny&1) * packetbits);
      return grid_demux(nbits,nx,(ny+1)/2,dx,dy*2,demux.ci,abits+1) + demux;
    }
  }
  

  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr auto grid_mux_preselect(u64 nbits, u64 nx, u64 ny, f64 dx/*um*/, f64 dy/*um*/, bool init=1, bool pol=0)
  {
    // each input data (nbits bits) is tagged with a "select" bit, which has been precomputed
    // all inputs but one have the tag bit reset
    // the tag propagates along with the data and is used as select signal for tristate inverters
    // the output tag is the ORing of the two input tags
    // returns two circuits: one for the tag, one for the data
    assert(nbits!=0);
    assert(nx!=0 && ny!=0);
    assert(dx>=0 && dy>=0);
    if (nx==1 && ny==1) {
      circuit data = (pol)? (inv{}.make(INVCAP) * nbits) : circuit{};
      return std::array{circuit{},data};
    }
    f64 psel = 1./(nx*ny); // a single select bit is set
    f64 bias = proba_bias(psel*0.5); // pseudo bias
    mux_inv_tri muxi{2}; // MUX (two tristate inverters)
    inv_tri itri; // tristate inverter
    f64 selcap = muxi.icap<1>() * nbits; // FIXME: wire not modeled
    f64 cselcap = muxi.icap<2>() * nbits; // FIXME: wire not modeled
    circuit locsel = buffer<SE,SMAX>(selcap,false,1,bias) | buffer<SE,SMAX>(cselcap,true,1,bias);
    if (init) {
      init = false;
      auto rest = grid_mux_preselect<SE,SMAX>(nbits,nx,ny,dx,dy,init,pol^1);
      circuit data = itri.make(rest[1].ci) * (nx*ny*nbits) + rest[1];
      circuit sel = locsel * (nx*ny*nbits) | rest[0];
      return std::array{sel,data};
    } else if ((nx>=2 && dx<=dy) || ny==1) {
      // reduce horizontally
      auto rest = grid_mux_preselect<SE,SMAX>(nbits,(nx+1)/2,ny,dx*2,dy,init,pol^1);
      circuit mux2 = muxi.make(rest[1].ci,1,bias);
      circuit mux1 = itri.make(rest[1].ci,1,bias);
      circuit datawire = wire<SE,SMAX>(dx/2,false,muxi.icap<0>(),false,bias);
      circuit data = (datawire * nx + (mux2 * (nx/2) || mux1 * (nx&1))) * (ny*nbits);
      data = data + rest[1];
      circuit tagcomp = nor{2}.make(rest[0].ci,1,bias);
      auto selwire = [=](bool cpl) {return wire<SE,SMAX>(dx/2,cpl,locsel.ci+tagcomp.ci,false,bias);};
      circuit selwires = selwire(true) * (2*(nx/2)*ny) || selwire(false) * ((nx&1)*ny);
      circuit sel = selwires + ((tagcomp * ((nx/2)*ny) + rest[0]) | locsel * (nx*ny));
      return std::array{sel,data};
    } else {
      // reduce vertically
      assert(ny>=2);
      auto rest = grid_mux_preselect<SE,SMAX>(nbits,nx,(ny+1)/2,dx,dy*2,init,pol^1);
      circuit mux2 = muxi.make(rest[1].ci,1,bias);
      circuit mux1 = itri.make(rest[1].ci,1,bias);
      circuit datawire = wire<SE,SMAX>(dy/2,false,muxi.icap<0>(),false,bias);
      circuit data = (datawire * ny + (mux2 * (ny/2) || mux1 * (ny&1))) * (nx*nbits);
      data = data + rest[1];
      circuit tagcomp = nor{2}.make(rest[0].ci,1,bias);
      auto selwire = [=](bool cpl) {return wire<SE,SMAX>(dy/2,cpl,locsel.ci+tagcomp.ci,false,bias);};
      circuit selwires = selwire(true) * (2*nx*(ny/2)) || selwire(false) * ((ny&1)*nx);
      circuit sel = selwires + ((tagcomp * (nx*(ny/2)) + rest[0]) | locsel * (nx*ny));
      return std::array{sel,data};
    }
  }


  template<u64 N, circuit OP, u64 D>
  inline constexpr circuit parallel_prefix = []() {
    // OP = associative 2-input operation, each input is D bits, N = total inputs
    // TODO: wires
    // FIXME: switching activity depends on OP
    static_assert(D!=0 && N!=0);
    if constexpr (N<=1) {
      return circuit{};
    } else {
      // Ladner-Fischer tree
      constexpr u64 K = (std::has_single_bit(N))? N/2 : std::bit_floor(N);
      constexpr circuit right = parallel_prefix<K,OP,D> + buffer(OP.ci*(N-K),false) * D;
      constexpr circuit left = parallel_prefix<N-K,OP,D>;
      return (left || right) + OP * (N-K);
    }
  }();
  

  template<u64 N>
  inline constexpr circuit priority_encoder = []() {
    // FIXME: not optimized
    static_assert(N!=0);
    if constexpr (N==1) {
      return circuit{};
    } else {
      constexpr circuit first = inv{}.make(INVCAP) * (N-1);
      constexpr circuit op = nanding(2,INVCAP) || noring(2,INVCAP);
      constexpr circuit last = noring(2,INVCAP) * (N-1);
      return first + parallel_prefix<N,op,2> + last;
    }
  }();


  template<bool INCR=false>
  constexpr circuit half_adder(f64 co)
  {
    if constexpr (INCR) {
      return inv{}.make(co); // a+1
    } else {
      return xor2(co) | anding(2,co); // a+b
    }
  }

  
  template<bool INCR=false>
  constexpr circuit full_adder(f64 co)
  {
    if constexpr (INCR) {
      return half_adder(co); // a+1+carry
    } else {
      return parity(3,co) | majority(co); // a+b+carry
    }
  }

  
  template<u64 N, bool INCR=false, bool CARRYIN=false>
  constexpr circuit adder_ks(f64 co)
  {
    // Kogge-Stone adder (radix-2)
    // TODO: wire capacitance
    static_assert(N!=0);
    if constexpr (N==1) {
      return (CARRYIN)? full_adder<INCR>(co) : half_adder<INCR>(co);
    }
    constexpr u64 depth = std::bit_width(N-1);
    basic_gate G[2] = {or_nand{},and_nor{}}; // inverting generate gate
    basic_gate P[2] = {nor{2},nand{2}}; // inverting propagate gate
    circuit sum = xor2(co) * (N-(CARRYIN^1)); // final sum
    circuit bws = (INCR)? circuit{} : xor2(sum.ci) * N; // bitwise sum (xnor for even depth, same circuit)
    circuit bwg; // bitwise generate
    circuit bwp; // bitwise propagate
    if constexpr (! INCR) {
      if constexpr (N==2) {
	bwg = nand{2}.make(G[0].icap<1>()+sum.ci) * N;
	bwp = nor{2}.make(G[0].icap<1>()+sum.ci) * (N-1);
      } else {
	bwg = nand{2}.make(G[0].icap<0>()+G[0].icap<1>()) * N;
	bwp = nor{2}.make(G[0].icap<1>()+2*P[0].icap()) * (N-1);
      }
    }
    circuit carryout;
    if (depth & 1) {
      carryout.ci = co;
    } else {
      carryout = inv{}.make(co);
    }
    circuit tree = bwp | bwg;
    for (int i=0; i<depth; i++) {
      u64 nd = 1<<i; // generate bits already calculated
      u64 ng = N-nd; // G cells
      u64 np = N-min(N,2*nd); // P cells
      u64 n2 = N-min(N,3*nd); // cells with 2 consumers
      u64 n1 = ng - n2; // cells with single consumer
      u64 np1 = min(n1,np); // P cells with single consumer
      if (i < (depth-1)) {
	f64 cog1 = G[(i+1)&1].icap<0>();
	f64 cog2 = G[(i+1)&1].icap<0>() + G[(i+1)&1].icap<1>();
	f64 cop1 = G[(i+1)&1].icap<1>() + P[(i+1)&1].icap();
	f64 cop2 = G[(i+1)&1].icap<1>() + 2 * P[(i+1)&1].icap();;
	circuit g = (G[i&1].make(cog1) * n1) || (G[i&1].make(cog2) * (ng-n1));
	circuit p = (P[i&1].make(cop1) * np1) || (P[i&1].make(cop2) * (np-np1));
	tree = tree + (g | p);
      } else {
	// last generate stage
	assert(np==0 && np1==0 && ng==n1);
	f64 cog1 = (n1==1)? carryout.ci : sum.ci;
	tree = tree + G[i&1].make(cog1) * n1;
      }
    } 
    return (tree | bws) + (sum || carryout);
  }
  
  // ###########################


  template<u64 N>
  struct flipflops {
    static constexpr circuit oneflop = []() {
      // edge-triggered flip-flop (2 latches)
      // clocking modeled separately
      inv i;
      inv_tri t;
      f64 co = INVCAP; // we do not care about the delay anyway
      f64 cfeedback1 = i.icap() + t.cp;
      f64 cfeedback2 = 2*i.icap() + t.cp;
      circuit tri1 = t.make(cfeedback1);
      circuit tri2 = t.make(cfeedback2);
      circuit forward = tri1 + i.make(2*t.icap()) + tri2 + i.make(co);
      circuit backward = tri1 || (i.make(t.icap()) + tri2);;
      return forward || backward;
    } ();

    static constexpr f64 height = SRAM_CELL.bitline_length * oneflop.t / SRAM_CELL.transistors; // um
    static constexpr f64 width = SRAM_CELL.wordline_length * N; // um

    static constexpr circuit flops = oneflop * N;
    
    static constexpr circuit clocking = []() {
      inv_tri t;
      f64 phicap1 = t.icap<1>() * 4;
      f64 phicap2 = t.icap<2>() * 4;
      circuit clock1 = buffer<4.>(phicap1,true); // complementary clock, local to each flop
      circuit clock2 = wire<4.>(width,false,(phicap2+clock1.ci)*N);
      return clock2 + clock1 * N;
    } ();

    static constexpr u64 xtors = flops.t + clocking.t;
    
    // clocking.e corresponds to a transition probability of 0.5
    // actually, the clock transitions with probability 1 twice per write
    static constexpr f64 write_energy_fJ = flops.e + clocking.e * 4;
  };

  
  // ###########################

  template<typename T>
  struct sram_common {
    static constexpr u64 NBITS = T::num_bits();
    static constexpr u64 XTORS = T::num_transistors();
    static constexpr f64 LATENCY = T::read_latency(); // ps
    static constexpr f64 EREAD = T::read_energy(); // fJ
    static constexpr f64 EWRITE = T::write_energy(); // fJ
    static constexpr f64 WIDTH = T::array_width(); // um
    static constexpr f64 HEIGHT = T::array_height(); // um
    // the following cost function is arbitrary (prioritizes latency over energy and reads over writes)
    static constexpr f64 COST = (10*EREAD+EWRITE) * pow(LATENCY,3);

    static void print(std::string s = "", std::ostream & os = std::cout)
    {
      os << s;
      os << std::setprecision(4);
      os << "bits: " << NBITS;
      os << " ; xtors: " << XTORS;
      os << " ; ps: " << LATENCY;
      os << " ;  fJ R|W: " << EREAD << "|" << EWRITE;
      os << " ;  um W|H: " << WIDTH << "|" << HEIGHT;
      //os << " ; area (um^2): " << WIDTH*HEIGHT;
      os << std::endl;
    }
  };
  

  struct sram_null {
    static constexpr u64 num_bits() {return 0;}
    static constexpr u64 num_transistors() {return 0;}
    static constexpr f64 read_latency() { return 0;}
    static constexpr f64 read_energy() {return 0;}
    static constexpr f64 write_energy() {return 0;}
    static constexpr f64 array_width() {return 0;}
    static constexpr f64 array_height() {return 0;}
  };


  template<u64 N, u64 M, u64 D>
  struct sram_bank {};

  
  template<u64 N, u64 M>
  struct sram_bank<N,M,0> : sram_common<sram_null>
  {
    static constexpr circuit ABUS {};
    static constexpr circuit WBUS {};
  };

  
  template<u64 N, u64 M, u64 D> requires (D<=M)
  struct sram_bank<N,M,D> : sram_common<sram_bank<N,M,D>> {
    // single R/W port, N wordlines, M cells per wordline
    // D = data width, M multiple of D, M/D power of 2
    static_assert(N!=0 && D!=0);
    static_assert((M%D)==0);
    static_assert(std::has_single_bit(M/D));
    
    // TODO: precharge circuit
    // TODO: double wordline, flying bitline (Chang, IEEE JSSC, jan 2021)
    // assume single-fin transistors (Chang, IEEE JSSC, jan 2021) ==> need write assist (TODO)
    
    static constexpr f64 WLCAP_fF = SRAM_CELL.wordline_capacitance;
    static constexpr f64 WLCAP = WLCAP_fF / CGATE_fF;
    static constexpr f64 WLCAP_pF = WLCAP * CGATE_pF;
    static constexpr f64 BLCAP_fF = SRAM_CELL.bitline_capacitance;
    static constexpr f64 BLCAP = BLCAP_fF / CGATE_fF;
    static constexpr f64 BLCAP_pF = BLCAP * CGATE_pF;
    static constexpr f64 WLRES = SRAM_CELL.wordline_resistance; // Ohm
    static constexpr f64 BLRES = SRAM_CELL.bitline_resistance; // Ohm

    // sense amplifier (SA) = latch type (cross coupled inverters)
    // SAVBLMIN value taken from Amrutur & Horowitz, IEEE JSSC, feb. 2000
    // TODO: SACAPMAX is currently a random number
    static constexpr f64 SACAPMAX = 30*INVCAP; // maximum SA input capacitance relative to CGATE
    static constexpr f64 SACAPMIN = INVCAP * (1+PINV); // minimum SA input capacitance relative to CGATE
    static_assert(SACAPMIN<=SACAPMAX);
    static constexpr f64 SAVBLMIN = 0.1; // minimum bitline voltage swing (V)
    // SA intrinsic delay
    static constexpr f64 SADELAY = 2 * inv{}.make(4*INVCAP).d; // 2 FO4 (Amrutur & Horowitz)
    static constexpr f64 XCSA = 0.4; // SA input capacitance relative to bitline capacitance
    
    // assumptions:
    //   * SA offset voltage follows Pelgrom's square-root law (Pileggi, CICC 2008)
    //   * bitline swing proportional to SA offset voltage (Abu-Rahma, CICC 2011)
    // ==> bitline swing inversely proportional to sqrt of SA input capacitance Csa (Kim, IEEE JSSC april 2023)
    static constexpr f64 SACAP = max(SACAPMIN,min(SACAPMAX,N*BLCAP*XCSA)); // sense amp capacitance relative to CGATE
    static constexpr f64 SACAP_pF = SACAP * CGATE_pF;
    static constexpr f64 SASCALE = SACAP/SACAPMIN;
    static_assert(SASCALE>=1);
    
    static constexpr f64 BLSWING = SAVBLMIN * sqrt(SACAPMAX/SACAP); // bitline voltage swing (V)

    // drive wordline from the mid point (Amrutur & Horowitz)
    static constexpr f64 WLRC = wire_res_delay((M/2)*WLRES, (M/2)*WLCAP_pF);

    static constexpr f64 BLRC = []() { // bitline RC delay (ps)
      // https://files.inria.fr/pacap/michaud/rc_delay.pdf
      // V = Vdd - (Idsat/2) * (t-t0) / (Cbl+Csa)
      // t0 = Rbl*Cbl*(Cbl+3Csa)/(6Cbl+6Csa)
      f64 CBL = N * BLCAP_pF; // fF
      f64 RBL = N * BLRES; // Ohm
      f64 T0 = RBL * CBL * (CBL+3*SACAP_pF) / (6*(CBL+SACAP_pF)); // ps
      return T0;     
    } ();
    
    static constexpr f64 BLDELAY = []() { // bitline total delay (ps)
      // model cell drive as current source (Amrutur & Horowitz)
      // neglect resistance of sense-amp isolation transistor
      // assume current is approximately Idsat/2 (2 nMOS in series)
      f64 CBL = N * BLCAP_pF; // fF
      return BLRC + (CBL+SACAP_pF) * BLSWING / (IDSAT/2);
    } ();

    // gates layout in peripheric logic is constrained by the wordline/bitline pitch
    // not sure to what extent this limits the gate size (TODO?)
    static constexpr u64 SMAX = 40;
    static constexpr f64 SEFF = 12.;
    
    static constexpr circuit RDEC = decode2<SEFF,SMAX>(N,M*WLCAP,SRAM_CELL.bitline_length); // row decoder

    // TODO: SA inverters are skewed
    // TODO: SA footer transistor (big capacitance)
    // one SA per column ==> assume SA does not impact column pitch (TODO?)
    static constexpr circuit SA = inv{}.make(SACAP,SASCALE) * 2 * M; // sense amplifiers
    
    // column read MUX (if D<M) after SA (wire capacitances not modeled, TODO)
    static constexpr auto CMUX = (D<M)? mux<SEFF,SMAX>(M/D,D,INVCAP) : std::array<circuit,2>{};

    static constexpr auto WDR = []() { // write driver
      std::array<circuit,2> wdr;
      circuit BLBUF = buffer<SEFF,SMAX>(N*BLCAP,false) | buffer<SEFF,SMAX>(N*BLCAP,true);
      // assume last inverter of buffer is tristate (impact on delay not modeled, TODO)
      if constexpr (D<M) {
	// data bits are interleaved
	// separate column decoder drives tristates
	constexpr u64 DEMUX = M/D;
	assert(DEMUX >= 2);
	f64 CTRI = 2 * N*BLCAP * TAU_ps / BLBUF.d; // input cap of tristate select (FIXME)
	if constexpr (DEMUX <= 16) {
	  // decode outputs run parallel to wordlines
	  f64 CSEL = CTRI*D + WLCAP*M;
	  circuit BUFTRI = buffer<SEFF>(CSEL,false) | buffer<SEFF>(CSEL,true);
	  circuit CDEC = decode2<SEFF>(DEMUX,BUFTRI.ci);
	  wdr[0] = CDEC + BUFTRI * DEMUX;
	  wdr[0].e += energy_fJ(WLCAP*M*CGATE_fF,VDD) * 2; // two wires switch (worst case)
	} else {
	  // predecode outputs run parallel to wordlines, AND2 gates are replicated (D replicas)
	  circuit BUFTRI = buffer<SEFF>(CTRI,false) | buffer<SEFF>(CTRI,true);
	  circuit CDEC = decode2_rep<SEFF>(DEMUX, BUFTRI.ci, SRAM_CELL.wordline_length, D);
	  wdr[0] = CDEC + BUFTRI * M;
	}
	wdr[1] = buffer<SEFF>((BLBUF.ci+WLCAP)*DEMUX,false);
	wdr[1].e += energy_fJ(WLCAP*DEMUX*CGATE_fF,VDD) * 0.5/*switch proba*/;
      }
      wdr[1] = wdr[1] * D + BLBUF * M;
      return wdr;
    } ();

    static constexpr f64 COLSELCAP = CMUX[0].ci + WDR[0].ci;
    static constexpr circuit ABUS = (buffer<SEFF,SMAX>(RDEC.ci,false) * std::bit_width(N-1)) || (buffer<SEFF,SMAX>(COLSELCAP,false) * std::bit_width(M/D-1));   
    static constexpr const circuit &WBUS = WDR[1];
    static constexpr circuit WLSEL = ABUS + RDEC;
    
    static constexpr f64 EWCL = D * 2 * inv{}.make(INVCAP).e; // cell energy per write
    static constexpr f64 EWL = 2 * energy_fJ(WLCAP_fF*M,VDD); // one wordline switches on and off
    // currently, assume perfect bitline voltage clamping (TODO?)
    // half-selected bitlines (D<M) consume energy on read (neglect on write, TODO?)
    static constexpr f64 EBLR = M * energy_fJ(N*BLCAP_fF,BLSWING) * 2; // bitline read + precharge
    static constexpr f64 EBLW = D * energy_fJ(N*BLCAP_fF,VDD) * 2; // bitline write + precharge
    
    static constexpr u64 CELL_XTORS = SRAM_CELL.transistors * N * M;
    static constexpr u64 PERI_XTORS = WLSEL.t + SA.t + WDR[0].t + WDR[1].t + CMUX[0].t + CMUX[1].t;
    static constexpr f64 PERI_GATECAP = WLSEL.cg + SA.cg + WDR[0].cg + WDR[1].cg + CMUX[0].cg + CMUX[1].cg;
    static constexpr f64 CELLS_AREA = SRAM_CELL.area * N * M; // um^2
    static constexpr f64 PERI_AREA = PERI_GATECAP / RANDOM_LOGIC_DENSITY; // um^2
    static constexpr f64 AREA = CELLS_AREA + PERI_AREA; // um^2
    static constexpr f64 CELLS_WIDTH = SRAM_CELL.wordline_length * M; // um
    static constexpr f64 CELLS_HEIGHT = SRAM_CELL.bitline_length * N; // um
    // assume that the overall aspect ratio equals the cells aspect ratio
    static constexpr f64 ASPECT_RATIO = CELLS_WIDTH / CELLS_HEIGHT;

    static constexpr u64 num_bits() {return N * M;}
    static constexpr u64 num_transistors() {return CELL_XTORS + PERI_XTORS;}
    static constexpr f64 array_width() {return sqrt(AREA*ASPECT_RATIO);} // um
    static constexpr f64 array_height() {return sqrt(AREA/ASPECT_RATIO);} // um

    static constexpr f64 read_latency()
    {
      return max(WLSEL.d + WLRC + BLDELAY + SADELAY, CMUX[0].d) + CMUX[1].d; // ps
    }

    static constexpr f64 read_energy()
    {
      return WLSEL.e + EWL + EBLR + SA.e + CMUX[0].e + CMUX[1].e; // fJ
    }

    static constexpr f64 write_energy()
    {
      return WLSEL.e + EWL + WDR[0].e + WDR[1].e + EBLW + EWCL; // fJ
    }

    static void print2(std::string s = "", std::ostream & os = std::cout)
    {
      os << s;
      WLSEL.print("wordline select: ",os);
      //SA.print("sense amps: ",os);
      if constexpr (D<M) {
	CMUX[0].print("column mux select: ",os);
	CMUX[1].print("column mux data: ",os);
      }
      WDR[0].print("write driver select: ",os);
      WDR[1].print("write driver data: ",os);
      os << std::setprecision(4);
      os << "bitline read voltage swing (V): " << BLSWING << std::endl;
      os << "wordline RC delay (ps): " << WLRC << std::endl;
      os << "bitline RC delay (ps): " << BLRC << std::endl;
      os << "bitline total delay (ps): " << BLDELAY << std::endl;
      os << "wordline energy (fJ): " << EWL << std::endl;
      os << "bitlines energy R/W (fJ): " << EBLR << " / " << EBLW << std::endl;
      os << "cells write energy: " << EWCL << std::endl;
      os << "sense amp scale: " << SASCALE << std::endl;
      os << "sense amp delay (ps): " << SADELAY << std::endl;
      os << "sense amps energy (fJ): " << SA.e << std::endl;
      os << "cell transistors: " << CELL_XTORS << std::endl;
      os << "periphery transistors: " << PERI_XTORS << std::endl;
      //os << "width (um): " << sram_bank::WIDTH << std::endl;
      //os << "height (um): " << sram_bank::HEIGHT << std::endl;
      os << "array efficiency: " << CELLS_AREA/AREA << std::endl;
    }
    
  };
  

  template<u64 N, u64 M, u64 D> requires (D>M)
  struct sram_bank<N,M,D> : sram_common<sram_bank<N,M,D>> {
    // single R/W port, N wordlines, M cells per wordline
    // split data over several banks, access all banks in parallel
    static_assert(N!=0 && M!=0);
    static constexpr u64 NB = D/M; // number of M-wide banks
    static constexpr u64 R = D%M; // extra bank provides R remaining bits
    using BANK = sram_bank<N,M,M>;
    using BANKR = sram_bank<N,R,R>;
    
    // send address to all banks
    static constexpr u64 ADDRESS_BITS = std::bit_width(N-1);
    static constexpr f64 LENGTH = BANK::WIDTH * NB + BANKR::WIDTH - 0.5 * (BANK::WIDTH + BANKR::WIDTH);
    static constexpr f64 LOADCAP = BANK::ABUS.ci * NB + BANKR::ABUS.ci;
    static constexpr circuit AWIRE = wire(LENGTH/2,false,LOADCAP,true);
    static constexpr circuit ABUS = (AWIRE | AWIRE) * ADDRESS_BITS; // address bus
    static constexpr circuit WBUS = BANK::WBUS * NB || BANKR::WBUS;
    
    static constexpr u64 num_bits() {return BANK::NBITS * NB + BANKR::NBITS;}
    static constexpr u64 num_transistors() {return BANK::XTORS * NB + BANKR::XTORS + ABUS.t;}
    static constexpr f64 read_latency() { return ABUS.d + BANK::LATENCY;}
    static constexpr f64 read_energy() {return ABUS.e + BANK::EREAD * NB + BANKR::EREAD;}
    static constexpr f64 write_energy() {return ABUS.e + BANK::EWRITE * NB + BANKR::EWRITE;}
    static constexpr f64 array_width() {return BANK::WIDTH * NB + BANKR::WIDTH;}
    static constexpr f64 array_height() {return BANK::HEIGHT;}

    static void print2(std::string s = "", std::ostream & os = std::cout)
    {
      os << s;
      ABUS.print("ABUS: ",os);
      BANK::print("BANK: ",os);
      if constexpr (R!=0) {
	BANKR::print("BANKR: ",os);
      }
    }
  };


  template<u64 BX, u64 BY, u64 N, u64 M, u64 D>
  struct sram_array : sram_common<sram_array<BX,BY,N,M,D>> {
    // BX x BY banks
    // Accessed data (D bits) lies in single bank
    // Bank is selected with XY-decoder before being accessed
    // Assume non-selected banks do not consume energy (TODO?)
    // Assume unidirectional interconnects (separate read/write data wires... FIXME: is this OK?)
    // Area of inter-bank wiring and logic is not modeled (TODO)
    static_assert(N!=0 && M!=0 && D!=0);
    using BANK = sram_bank<N,M,D>;
    static constexpr u64 NB = BX * BY; // number of banks
    static_assert(NB!=0);
    // N must be a power of 2 unless there is a single bank
    static_assert(NB==1 || std::has_single_bit(N));
    // BY must be a power of 2 unless BX=1
    static_assert(BX==1 || std::has_single_bit(BY));
    static constexpr f64 HW = BANK::WIDTH * ((BX-1)/2); // half width
    static constexpr f64 HH = BANK::HEIGHT * ((BY-1)/2); // half height
    static constexpr circuit SELG = (BX==1 || BY==1)? circuit{} : anding(2,INVCAP,1,(1./BX+1./BY)/2);
    static constexpr circuit XSEL = wire(HW,false,SELG.ci*(BY/2),true,1./BX) | wire(HW,false,SELG.ci*(BY-BY/2),true,1./BX);
    static constexpr circuit YSEL = wire(HH,false,SELG.ci*(BX/2),true,1./BY) | wire(HH,false,SELG.ci*(BX-BX/2),true,1./BY);
    static constexpr circuit XDEC = (BX==1)? circuit{} : decode2(BX,YSEL.ci,BANK::WIDTH) + YSEL * BX;
    static constexpr circuit YDEC = (BY==1)? circuit{} : decode2(BY,XSEL.ci,BANK::HEIGHT) + XSEL * BY;
    static constexpr circuit SEL = (XDEC || YDEC) + SELG * (BX*BY); // bank select
    static constexpr u64 ABITS = (N>=2)? std::bit_width(N-1) : 1; // local (bank) address bits
    static constexpr circuit ATREE = grid_demux(ABITS,BX,BY,BANK::WIDTH,BANK::HEIGHT,BANK::ABUS.ci); // address tree
    static constexpr circuit WTREE = grid_demux(D,BX,BY,BANK::WIDTH,BANK::HEIGHT,BANK::WBUS.ci); // data write tree
    static constexpr auto RTREE = grid_mux_preselect(D,BX,BY,BANK::WIDTH,BANK::HEIGHT); // read TREE
    static constexpr circuit ACC = SEL || ATREE;
    static constexpr circuit READ = []() {
      circuit bankread;
      bankread.d = BANK::LATENCY;
      bankread.e = BANK::EREAD;
      // bank access starts after bank select signal has been broadcast
      return ACC + (RTREE[0] || (bankread + RTREE[1]));
    }();
    
    static void print2(std::string s = "", std::ostream & os = std::cout)
    {
      ATREE.print("ATREE: ",os);
      SEL.print("SEL: ",os);
      RTREE[0].print("RTREE select: ",os);
      RTREE[1].print("RTREE data: ",os);
      WTREE.print("WTREE: ",os);
      BANK::print("NODE: ",os);
      //BANK::print2("",os);
    }
    
    static constexpr u64 num_bits() {return NB * BANK::NBITS;}
    static constexpr u64 num_transistors() {return NB * BANK::XTORS + ACC.t + WTREE.t + RTREE[0].t + RTREE[1].t;}
    static constexpr f64 read_latency() {return READ.d;}
    static constexpr f64 read_energy() {return READ.e;}
    static constexpr f64 write_energy() {return ACC.e + WTREE.e + BANK::EWRITE;}
    static constexpr f64 array_width() {return BANK::WIDTH * BX;}
    static constexpr f64 array_height() {return BANK::HEIGHT * BY;}
  };


  constexpr bool ok_config(u64 E/*entries*/, u64 D/*data bits*/, u64 MAXN, u64 MAXM)
  {
    bool ok = E!=0 && D!=0; // mandatory
    ok = ok && std::has_single_bit(MAXN); // mandatory: MAXN must be a power of 2
    // prune off configs that are likely bad
    u64 banksize = MAXN * max(D,MAXM);
    ok = ok && (16*E*D <= banksize*banksize || MAXN>=1024 && MAXM>=512);
    return ok;
  }
  

  template<u64 E, u64 D, u64 MAXN, u64 MAXM>
  struct sram_banked;

  
  template<u64 E, u64 D, u64 MAXN, u64 MAXM> requires (ok_config(E,D,MAXN,MAXM))
    struct sram_banked<E,D,MAXN,MAXM> : sram_common<sram_banked<E,D,MAXN,MAXM>> {
    static_assert(E!=0 && D!=0 && std::has_single_bit(MAXN));
    static constexpr bool ok = true;

    static constexpr auto banking = []<u64 R/*rows*/, u64 W/*row width*/>() {
      constexpr u64 N = MAXN;
      static_assert(R>=N);
      static_assert(W>=MAXM/2);
      constexpr f64 BIAR = f64(N) / (W*SRAM_CELL.aspect_ratio); // bank inverse aspect ratio
      constexpr f64 FACTOR = max(2,BIAR*BIAR);
      static_assert(FACTOR>=1);
      if constexpr (R <= W * SRAM_CELL.aspect_ratio * FACTOR) {
	// single column
	return std::tuple{1/*BX*/,(R+N-1)/N/*BY*/};
      } else {
	// multiple columns (BY=2^k)
	// we want squarish shape
	constexpr f64 S = sqrt(R*W*SRAM_CELL.aspect_ratio); // square side
	static_assert(S>=N);
	static_assert(S<=R);
	constexpr u64 BY = min(to_pow2(S/N),std::bit_floor(R/N));
	static_assert(BY>=1 && std::has_single_bit(BY));
	constexpr u64 BX = (R+BY*N-1)/(BY*N);
	static_assert(BX>=1);
	return std::tuple{BX,BY};
      }
    };
    
    static constexpr auto params = []() {
      if constexpr (D > MAXM) {
	// data width = bank width
	constexpr u64 M = MAXM;
	if constexpr (E <= MAXN) {
	  // single bank
	  return std::tuple{1/*BX*/,1/*BY*/,E/*N*/,M};
	} else {
	  // multiple banks (N=MAXN)
	  auto [BX,BY] = banking.template operator()<E,D>();
	  return std::tuple{BX,BY,MAXN,M};
	}
      } else {
	// D <= MAXM
	constexpr u64 MAXK = std::bit_floor(MAXM/D);
	if constexpr (E <= MAXN*MAXK) {
	  // single bank
	  u64 K = MAXK;
	  u64 N = (E+K-1)/K;
	  // we want squarish shape
	  while (K*D*SRAM_CELL.aspect_ratio > 2*N && K>=2 && (E+K/2-1)/(K/2)<=MAXN) {
	    K /= 2;
	    N = (E+K-1)/K;
	  }
	  u64 M = K * D;
	  return std::tuple{1/*BX*/,1/*BY*/,N,M};
	} else {
	  // multiple banks (N=MAXN)
	  constexpr u64 M = MAXK * D;
	  constexpr u64 EK = (E+MAXK-1) / MAXK;
	  auto [BX,BY] = banking.template operator()<EK,M>();
	  return std::tuple{BX,BY,MAXN,M};
	}
      }
    } ();

    static constexpr u64 BX = std::get<0>(params);
    static constexpr u64 BY = std::get<1>(params);
    static constexpr u64 N = std::get<2>(params);
    static constexpr u64 M = std::get<3>(params);
    
    using ARR = sram_array<BX,BY,N,M,D>;
    static_assert(ARR::NBITS >= E*D);

    static void print2(std::string s = "", std::ostream & os = std::cout)
    {
      os << s;
      os << "BX=" << BX;
      os << " BY=" << BY;
      os << " N=" << N;
      os << " M=" << M;
      os << std::endl;
      //ARR::print2("",os);
    }
    
    static constexpr u64 num_bits() {return ARR::NBITS;}
    static constexpr u64 num_transistors() {return ARR::XTORS;}
    static constexpr f64 read_latency() {return ARR::LATENCY;}
    static constexpr f64 read_energy() {return ARR::EREAD;}
    static constexpr f64 write_energy() {return ARR::EWRITE;}
    static constexpr f64 array_width() {return ARR::WIDTH;}
    static constexpr f64 array_height() {return ARR::HEIGHT;}
  };


  template<u64 E, u64 D, u64 MAXN, u64 MAXM> requires (! ok_config(E,D,MAXN,MAXM))
  struct sram_banked<E,D,MAXN,MAXM> : sram_common<sram_null> {
    static constexpr bool ok = false; 
  };

  
  template<typename T>
  concept sram_type = requires(T& x) {[]<u64 E, u64 D, u64 N, u64 M>(sram_banked<E,D,N,M>&){}(x);};

  template<sram_type T1, sram_type ...T>
  struct best_config;

  template<sram_type T>
  struct best_config<T> {
    using type = T;
  };

  template<sram_type T1, sram_type T2, sram_type ...T>
  struct best_config<T1,T2,T...> {
    using best = best_config<T2,T...>::type;
    static constexpr bool better = ! best::ok || (T1::ok && (T1::COST <= best::COST));
    using type = std::conditional_t<better,T1,best>;
  };

  template<u64 E, u64 D, u64 N, u64 ...M>
  using sram_bestM = best_config<sram_banked<E,D,N,M>...>::type;

  template<u64 E, u64 D, u64 ...N>
  using sram_bestN = best_config<sram_bestM<E,D,N,64,128,256,512>...>::type;
  
  template<u64 E, u64 D>
  using sram =sram_bestN<E,D,64,128,256,512,1024>;

  template<u64 E, u64 D0, u64 S, u64 ...I>
  void test_sram_D(std::integer_sequence<u64,I...> seq)
  {
    (sram<E,D0+I*S>::print("D="+std::to_string(D0+I*S)+": "),...);
  }

  template<u64 D, u64 E0, f64 X, u64 ...I>
  void test_sram_E(std::integer_sequence<u64,I...> seq)
  {
    (sram<llround(E0*pow(X,I)),D>::print("E="+std::to_string(llround(E0*pow(X,I)))+": "),...);
  }  

  
  // ###########################
  
  inline constexpr f64 OUTCAP = INVCAP;
  
  inline constexpr circuit INV = inv{}.make(OUTCAP);
  
  template<u64 N>
  inline constexpr circuit AND = anding(N,OUTCAP);
  
  template<u64 N>
  inline constexpr circuit OR = oring(N,OUTCAP);

  template<u64 N>
  inline constexpr circuit NAND = nanding(N,OUTCAP);
  
  template<u64 N>
  inline constexpr circuit NOR = noring(N,OUTCAP);

  template<u64 N> requires (N>=2)
  inline constexpr circuit XOR = parity(N,OUTCAP);
  
  template<u64 N> requires (N>=2)
  inline constexpr circuit XNOR = XOR<N>;
  
  template<u64 WIDTH>
  inline constexpr circuit ADD = adder_ks<WIDTH,false,false>(OUTCAP);

  template<u64 WIDTH>
  inline constexpr circuit INC = adder_ks<WIDTH,true,true>(OUTCAP);

  template<u64 WIDTH>
  inline constexpr circuit SUB = adder_ks<WIDTH,false,true>(OUTCAP);
  
  template<u64 N>
  inline constexpr circuit EQUAL = xnor2(AND<N>.ci) * N + AND<N>;

  template<u64 N>
  inline constexpr circuit NEQ = xor2(OR<N>.ci) * N + OR<N>;

  template<u64 N, u64 WIDTH> requires (N>=2)
  inline constexpr auto MUX = mux(N,WIDTH,OUTCAP);
  
  
  // ###########################

  template<typename T>
  struct base_impl {};

  template<arith T>
  struct base_impl<T> {using type = T;};
  
  template<valtype T>
  struct base_impl<T> {using type = T::type;};

  template<typename T>
  using base = base_impl<T>::type;

  template<typename T>
  concept ival = std::integral<T> || std::integral<base<T>>; 

  template<typename T>
  concept fval = std::floating_point<T> || std::floating_point<base<T>>; 
  
  template<typename T>
  inline constexpr u64 length = 0;

  template<valtype T>
  inline constexpr u64 length<T> = T::size;

  template<typename T>
  struct toval_impl {};
  
  template<arith T>
  struct toval_impl<T> {
    using type = val<sizeof(T)*8,T>;
  };

  template<valtype T>
  struct toval_impl<T> {
    using type = val<T::size,typename T::type>;
  };
  
  template<typename T>
  using toval = toval_impl<T>::type;
  
  template<typename T1, typename ...T>
  struct valt_impl {};

  template<typename T>
  struct valt_impl<T> {
    using type = toval<T>;
  };
  
  template<ival T1, ival T2>
  struct valt_impl<T1,T2> {
    using type = std::conditional_t<(length<T2> > length<T1>),toval<T2>,toval<T1>>;
  };

  template<fval T1, fval T2>
  struct valt_impl<T1,T2> {
    using type = std::conditional_t<(length<T2> > length<T1>),toval<T2>,toval<T1>>;
  };

  template<ival T1, fval T2>
  struct valt_impl<T1,T2> {
    using type = toval<T2>;
  };  

  template<fval T1, ival T2>
  struct valt_impl<T1,T2> {
    using type = toval<T1>;
  };

  template<typename T1, typename T2, typename ...T>
  struct valt_impl<T1,T2,T...> {
    using type = typename valt_impl<T1,typename valt_impl<T2,T...>::type>::type;
  };
  
  template<typename ...T>
  using valt = valt_impl<T...>::type;

  
  // ###########################

  class exec_control {
    template<u64,arith> friend class val;
    template<memdatatype,u64> friend class ram;
    friend class globals;
    template<u64 N, action T> friend void execute(val<N,u64>,T);
  private:
    bool active = true;
    u64 time = 0;
    exec_control(const exec_control &s) = default;
    exec_control& operator=(const exec_control &s) = default;
    void set_state(bool cond, u64 t=0) {active = cond; time = t;}
  public:
    exec_control() : active(true), time(0) {}
  } exec;

  
  // ###########################
  
  template<arith T>
  class global {
    friend class globals;
    friend class ::simulator;
    template<u64,arith> friend class reg;
    template<memdatatype,u64> friend class ram;
  private:
    T data = 0;
    global() : data(0) {}
    global(const global&) = default;
    global& operator= (const global& x) = default;
    void operator+= (T i) {data+=i;}
    void operator++(int) {data++;}
    operator T() const {return data;}
    void operator&() = delete;
  public:
    global(T x) : data(x) {}
    
    void print(std::string s = "", std::ostream & os = std::cout) const
    {
      os << s << data << std::endl;
    }
  };


  class globals {
    template<u64,arith> friend class val;
    template<u64,arith> friend class reg;
    template<valtype,u64> friend class arr;
    template<memdatatype,u64> friend class ram;
    friend class proxy;
    friend class ::simulator;
  public:
    global<u64> clock_cycle_ps;
    global<u64> storage;
    global<f64> energy_fJ;
    //global<u64> transistors;
  private:
    void update_storage(u64 n) {storage += n;}
    void update_energy(f64 e) {if (exec.active) energy_fJ += e;}

    void update_metrics(const circuit &c)
    {
      //transistors += c.t;
      update_energy(c.energy());
    }
  } panel;

  
  // ###########################

  
  template<u64 N, arith T = u64>
  class val {
    //static_assert(N!=0);
    static_assert(N <= sizeof(T)*8);
    static_assert(N==sizeof(T)*8 || std::integral<T>);

    template<u64,arith> friend class val;
    template<u64,arith> friend class reg;
    template<valtype, u64> friend class arr;
    template<memdatatype,u64> friend class ram;
    template<u64,u64,arith> friend class rom;
    friend class proxy;
    friend class ::simulator;
    
  private:
    T data = 0;
    u64 timing = 0; // time

    val() {};
    
    T fit(T x) requires std::integral<T>
    {
      return truncate<N>(x);
    }
    
    T get() const {return data;}

    T get() const requires std::signed_integral<T>
    {
      if constexpr (N < sizeof(T)*8) {
	assert((std::make_unsigned_t<T>(data) >> N) == 0);
      }
      T signbit = T(1) << (N-1);
      return (data ^ signbit) - signbit;
    }

    u64 time() const {return max(exec.time,timing);}

    void set_time(u64 t) {timing = t;}

    void operator= (const val& x)
    {
      if (! exec.active) return;
      data = x.get();
      timing = x.time();
    }

    void operator&() = delete;
    
  public:
    
    static constexpr u64 size = N;
    using type = T;

    static constexpr T maxval = []() {
      if constexpr (N == sizeof(T)*8) {
	return std::numeric_limits<T>::max();
      } else {
	static_assert(N < sizeof(T)*8);
	if constexpr (std::unsigned_integral<T>) {
	  return (T(1)<<N)-1;
	} else {
	  static_assert(std::signed_integral<T>);
	  return (T(1)<<(N-1))-1;
	}
      } 
    }();

    static constexpr T minval = []() {
      if constexpr (N == sizeof(T)*8) {
	return std::numeric_limits<T>::min();
      } else {
	static_assert(N < sizeof(T)*8);
	if constexpr (std::unsigned_integral<T>) {
	  return 0;
	} else {
	  static_assert(std::signed_integral<T>);
	  return -(T(1)<<(N-1));
	}
      } 
    }();    

    val(std::integral auto x, u64 t=0) requires std::integral<T> : data(fit(x))
    {
      set_time(t);
    }

    val(std::integral auto x, u64 t=0) : data(x)
    {
      set_time(t);
    }
    
    val(std::floating_point auto x, u64 t=0) : data(x)
    {
      static_assert(N == sizeof(T)*8);
      set_time(t);
    }

    val(const val& x) : val(x.get(),x.time()) {}
    
    val(const valtype auto& x) requires std::unsigned_integral<T> : val(to_unsigned(x.get()),x.time()) {}

    val(const valtype auto& x) : val(x.get(),x.time()) {}
    
    void print(bool t=true, std::string before="", std::string after="\n", std::ostream & os=std::cout) const
    {
      os << before << +get();
      if (t)
	os << " (t=" << time() << ")";
      os << after << std::flush;
    }
    
    void printb(bool t=true, std::string before="", std::string after="\n", std::ostream & os=std::cout) const
    {
      os << before;
      if constexpr (std::integral<T>) {
	os << std::bitset<N>(get());
      } else if constexpr (std::same_as<T,f32>) {
	os << std::bitset<N>(std::bit_cast<u32>(get()));
      } else if constexpr (std::same_as<T,f64>) {
	os << std::bitset<N>(std::bit_cast<u64>(get()));
      }
      if (t)
	os << " (t=" << time() << ")";
      os << after << std::flush;
    }

    val<1> operator[] (u64 i) requires std::unsigned_integral<T>
    {
      // no transistors
      // only for reading a bit, not for writing
      assert(i<N);
      return {(get() >> i) & 1, time()};
    }
    
    [[nodiscard]] val reverse() const requires std::unsigned_integral<T>
    {
      // no transistors
      return {reverse_bits(get()) >> (sizeof(T)*8-N), time()};
    }

    auto ones() const
    {
      // TODO: transistors & time
      auto n = std::popcount(truncate<N>(get()));
      return val<std::bit_width(N)> {n};
    }

    [[nodiscard]] val priority_encode() const requires std::unsigned_integral<T>
    {
      constexpr circuit c = priority_encoder<N>;
      u64 x = get();
      u64 y = x & (x^(x-1));
      panel.update_metrics(c);
      return {y,time()+c.delay()};
    }

    template<u64 M>
    auto make_array() const
    {
      // broadcast (TODO: wires)
      arr<valt<val>,M> a;
      for (u64 i=0; i<M; i++) {
	a[i] = *this;
      }
      return a;
    }
  };
  
  
  // ###########################

  class proxy {
  private:
    proxy() = delete;
    ~proxy() = delete;
    static auto get(const auto& x) {return x;}
    static auto get(const valtype auto& x) {return x.get();}
    static auto get(const arrtype auto& x) {return x.get();}
    static auto time(const auto& x) {return 0;}
    static auto time(const valtype auto& x) {return x.time();}
    static auto time(const arrtype auto& x) {return x.time();}

    static void set_time(valtype auto& x, arith auto t) {x.set_time(t);}
    static void update_metrics(const circuit &c) {panel.update_metrics(c);}
    
    template<valtype T1, valtype T2>
    static auto concatenate(const T1& x1, const T2& x2)
    {
      static_assert(std::unsigned_integral<typename T1::type>);
      static_assert(std::unsigned_integral<typename T2::type>);
      return val<x1.size+x2.size> {(proxy::get(x1) << x2.size) | proxy::get(x2)};
    }
    
    template<u64,arith> friend class reg;
    template<valtype,u64> friend class arr;
    template<memdatatype,u64> friend class ram;
    template<u64,u64,arith> friend class rom;
    
    template<valtype T1, valtype T2>
    friend val<1> operator== (const T1&, const T2&);

    template<valtype T1, arith T2>
    friend val<1> operator== (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend val<1> operator!= (const T1&, const T2&);

    template<valtype T1, arith T2>
    friend val<1> operator!= (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend val<1> operator> (const T1&, const T2&);

    template<valtype T1, arith T2>
    friend val<1> operator> (const T1&, const T2&);

    template<valtype T1, valtype T2>
    friend val<1> operator< (const T1&, const T2&);
    
    template<valtype T1, arith T2>
    friend val<1> operator< (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend val<1> operator>= (const T1&, const T2&);

    template<valtype T1, arith T2>
    friend val<1> operator>= (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend val<1> operator<= (const T1&, const T2&);

    template<valtype T1, arith T2>
    friend val<1> operator<= (const T1&, const T2&);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1,T2> operator+ (const T1&, const T2&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1,T2> operator+ (const T1&, const T2&);

    template<valtype T>
    friend valt<T> operator- (const T&);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1,T2> operator- (const T1&, const T2&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1,T2> operator- (const T1&, const T2&);
    
    template<std::integral T1, valtype T2> requires (ival<T2>)
    friend valt<T1,T2> operator- (const T1&, const T2&);

    template<valtype T1, valtype T2>  requires (ival<T1> && ival<T2>)
    friend valt<T1> operator<< (const T1&, const T2&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1> operator<< (const T1&, const T2&);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1> operator>> (const T1&, const T2&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1> operator>> (const T1&, const T2&);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1,T2> operator* (const T1&, const T2&);

    template<valtype T1, ival T2> requires (ival<T1>)
    friend valt<T1,T2> operator* (const T1&, const T2&);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1> operator/ (const T1&, const T2&);

    template<valtype T1, ival T2> requires (ival<T1>)
    friend valt<T1> operator/ (const T1&, const T2&);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T2> operator% (const T1&, const T2&);

    template<valtype T1, ival T2> requires (ival<T1>)
    friend valt<T2> operator% (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator& (const T1&, const T2&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator& (const T1&, const T2&);

    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator| (const T1&, const T2&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator| (const T1&, const T2&);
    
    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator^ (const T1&, const T2&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator^ (const T1&, const T2&);
    
    template<valtype T>
    friend valt<T> operator~ (const T&);

    template<valtype T1, valtype T2, valtype... T>
    friend auto concat(const T1&, const T2&, const T&...);

    template<valtype T>
    friend valt<T> select(val<1>,const T&,const T&);

    template<u64 N, action T>
    friend void execute(val<N,u64>,T);
  };
  
  // ###########################

  inline global<bool> storage_destroyed = false;
  
  
  template<u64 N, arith T = u64>
  class reg : public val<N,T> {
  public:
    using stg = flipflops<N>;
    
  private:
    void create()
    {
      assert(("all storage (reg,ram) must have the same lifetime",!storage_destroyed));
      panel.update_storage(N);
    }
    
  public:
    
    reg() : val<N,T>() {create();}

    reg(const std::convertible_to<T> auto &x) : val<N,T>(x) {create();}

    template<u64 M, arith U>
    reg(const val<M,U> &x) : val<N,T>(x) {create();}

    ~reg()
    {
      storage_destroyed = true;
    }

    void operator= (const auto &other)
    {
      val<N,T>::operator=(other);
      panel.update_energy(stg::write_energy_fJ);
    }

    reg(const reg &other)
    {
      *this = other;
      create();
    }

  };


  // ###########################

  template<typename T, typename X, typename Y>
  concept unaryfunc = requires (T f, X i) {{f(i)} -> std::convertible_to<Y>;};

  
  template<valtype T, u64 N>
  class arr {
    template<valtype,u64> friend class arr;
    template<u64,arith> friend class val;
    template<memdatatype,u64> friend class ram;
    friend class proxy;
  private:

    std::array<T,N> a {};

    void copy_from(const arr &other)
    {
      for (u64 i=0; i<N; i++) a[i] = other.a[i];
    }
    
    void operator= (const arr &other) requires (! regtype<T>)
    {
      copy_from(other);
    }

    void operator&() = delete;

    auto get() const
    {
      std::array<typename T::type,N> b;
      for (int i=0; i<N; i++) b[i] = a[i].get();
      return b;
    }
    
    u64 time() const
    {
      u64 t = 0;
      for (u64 i=0; i<N; i++) {
	t = max(t,a[i].time());
      }
      return t;
    }
    
    void set_time(u64 t)
    {
      for (u64 i=0; i<N; i++) a[i].set_time(t);
    }
    
  public:

    static constexpr u64 size = N;
    using type = T;

    arr() {}
    
    template<std::convertible_to<T> ...U>
    arr(U... args) : a{args...} {}
    
    arr(unaryfunc<u64,T> auto f)
    {
      for (u64 i=0; i<N; i++) {
	a[i] = f(i);
      }
    }

    template<std::convertible_to<T> U>
    arr(const std::array<U,N> &b)
    {
      for (u64 i=0; i<N; i++) {
	a[i] = b[i];
      }
    }

    arr(const arr &other) {copy_from(other);}    

    void operator= (const arr &other) requires (regtype<T>)
    {
      copy_from(other);
    }
    
    T get(u64 i) const
    {
      assert(i<N);
      return a[i];
    }
    
    T& operator [] (u64 i)
    {
      assert(i<N);
      return a[i];
    }

    void print(bool t=true, std::string before="", std::string after="\n", std::ostream & os=std::cout) const
    {
      for (u64 i=0; i<N; i++) a[i].print(t,before+std::to_string(i)+": ",after,os);
    }

    void printb(bool t=true, std::string before="", std::string after="\n", std::ostream & os=std::cout) const
    {
      for (u64 i=0; i<N; i++) a[i].printb(t,before+std::to_string(i)+": ",after,os);
    }

    auto concat() const requires std::unsigned_integral<typename T::type>
    {
      // element 0 is at rightmost position
      // no transistors (TODO: wires)
      static_assert(N!=0);
      u64 y = 0;
      for (i64 i=N-1; i>=0; i--) {
	y = (y<<T::size) | a[i].get();
      }
      return val<N*T::size> {y,time()};
    }

    [[nodiscard]] auto append(const std::convertible_to<valt<T>> auto& x) const
    {
      // no transistors
      arr<valt<T>,N+1> b;
      for (int i=0; i<N; i++) {
	b[i] = a[i];
      }
      b[N] = x;
      return b;
    }
    
    template<u64 M>
    val<M> extract(u64 pos) const requires std::unsigned_integral<typename T::type>
    {
      // no transistors
      static_assert(M!=0);
      if (pos >= T::size*N)
	return 0;
      int i = pos / T::size;
      int j = pos % T::size;
      u64 x = a[i].get() >> j;
      auto t = a[i].time();
      u64 nbits = T::size - j;
      i++;
      while (nbits<M && i<N) {
	x |= a[i].get() << nbits;
	t = max(t,a[i].time());
	nbits += T::size;
	i++;
      }
      return {x,t};
    }

    template<u64 W, u64 NBITS>
    auto make_array() const requires std::unsigned_integral<typename T::type>
    {
      // no transistors
      static_assert(W!=0 && NBITS!=0);
      constexpr u64 L = (NBITS+W-1)/W;
      static_assert(L!=0);
      arr<val<W>,L> b;
      for (int i=0; i<L; i++) {
	b[i] = extract<W>(i*W);
      }
      if constexpr (constexpr u64 R = NBITS % W; R != 0) {
	b[L-1] = b[L-1].get() % (1ull<<R);
      }
      return b;
    }   
    
    valt<T> xor_all() const
    {
      if constexpr (N>=2) {
	static constexpr circuit c = XOR<N> * T::size;
	panel.update_metrics(c);
	auto x = a[0].get();
	for (int i=1; i<N; i++) {
	  x ^= a[i].get();
	}
	u64 t = time() + c.delay();
	return {x,t};
      } else if constexpr (N==1) {
	// no transistors
	return a[0];
      } else {
	return 0;
      }
    }
  };


  // ###########################

  template<memdatatype T>
  struct rawdata {};

  template<memdatatype T> requires valtype<T>
  struct rawdata<T> {
    using type = T::type;
    static constexpr u64 width = T::size;
  };

  template<memdatatype T> requires arrtype<T>
  struct rawdata<T> {
    using type = std::array<typename T::type::type,T::size>;
    static constexpr u64 width = T::size * T::type::size;
  };  

  
  template<memdatatype T, u64 N>
  class ram {
    // bandwidth = one read per clock cycle
  public:
    using type = T;
    using valuetype = rawdata<T>::type;
    using static_ram = sram<N,rawdata<T>::width>;

  private:
    valuetype data[N];
    u64 next_read_time = 0;

    struct writeop {
      u64 addr{};
      valuetype dataval{};
      u64 t{};
      
      bool operator< (const writeop &rhs) {return t > rhs.t;}

      void commit(ram &mem) const
      {
	assert(addr<N);
	mem.data[addr] = dataval;
      }
    };
    
    std::vector<writeop> writes; // pending writes

    ram(const ram &) = delete;
    ram& operator= (const ram&) = delete;
    void operator& () = delete;

  public:
    
    ram()
    {
      panel.update_storage(static_ram::NBITS);
    }

    ~ram()
    {
      storage_destroyed = true;
      while (! writes.empty()) {
	writes[0].commit(*this);
	std::pop_heap(writes.begin(),writes.end());
	writes.pop_back();
      }
    }

    void write(const ival auto& address, const std::convertible_to<T> auto& dataval)
    {
      if (! exec.active) return;
      panel.update_energy(static_ram::EWRITE);
      u64 t = max(proxy::time(address),proxy::time(dataval)); // TODO: write latency?
      if (t <= next_read_time) {
	// write immediately
	u64 addr = proxy::get(address);
	assert(addr<N);
	data[addr] = proxy::get(dataval);
      } else {
	// write in the future
	writes.push_back({u64(proxy::get(address)),valuetype(proxy::get(dataval)),t});
	std::push_heap(writes.begin(),writes.end());
      }
    }

    T read(const ival auto& address)
    {
      if (! exec.active) return {};
      panel.update_energy(static_ram::EREAD);
      u64 t = max(proxy::time(address),next_read_time); // time at which the read starts
      while (! writes.empty() && writes[0].t <= t) {
	writes[0].commit(*this);
	std::pop_heap(writes.begin(),writes.end());
	writes.pop_back();
      }
      next_read_time = t + panel.clock_cycle_ps;
      u64 addr = proxy::get(address);
      assert(addr<N);
      t += std::llround(static_ram::LATENCY); // time at which the read completes
      T readval = data[addr];
      readval.set_time(t);
      return readval;
    }
    
    void reset()
    {
      if (! exec.active) return;
      panel.update_energy(static_ram::NBITS*2*INV.e);
      for (auto &v : data) v = {}; // instantaneous (FIXME?)
      writes.clear();
    }

    void print(std::string s = "", std::ostream & os = std::cout)
    {
      static_ram::print(s,os);
    }
  };

  
  // ###########################
  
  template<u64 N, u64 M, arith T = u64>
  class rom {
    val<M,T> a[N];
    
  public:

    template<std::convertible_to<T> U>
    rom(std::initializer_list<U> l)
    {
      assert(l.size() == N);
      int i = 0;
      for (auto e : l) a[i++] = e;       
    }

    template<std::convertible_to<T> U>
    rom(U (&x)[N])
    {
      for (int i=0; i<N; i++) a[i] = x[i];
    }

    val<M,T> read(const ival auto& address) const
    {
      auto i = proxy::get(address);
      assert(i<N);
      u64 t = 0;
      if constexpr (valtype<decltype(address)>) {
	t = max(t,address.time());
      }
      //TODO: energy and delay
      return {proxy::get(a[i]),t};
    }    
  };

    
  // ###########################

  // EQUALITY
  template<valtype T1, valtype T2>
  val<1> operator== (const T1& x1, const T2& x2)
  {
    static_assert(T1::size == T2::size);
    static constexpr circuit c = EQUAL<T1::size>;
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) == proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator== (const T1& x1, const T2& x2)
  {
    static constexpr circuit reduction = NOR<T1::size>;
    circuit c = INV * ones<T1::size>(x2) + reduction; // not constexpr
    proxy::update_metrics(c);
    auto t = proxy::time(x1) + c.delay();
    return {proxy::get(x1) == x2, t};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator== (const T1& x1, const T2& x2)
  {
    return x2 == x1;
  }

  // INEQUALITY
  template<valtype T1, valtype T2>
  val<1> operator!= (const T1& x1, const T2& x2)
  {
    static_assert(T1::size == T2::size);
    static constexpr circuit c = NEQ<T1::size>;
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) != proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator!= (const T1& x1, const T2& x2)
  {
    static constexpr circuit reduction = OR<T1::size>;
    circuit c = INV * ones<T1::size>(x2) + reduction; // not constexpr
    proxy::update_metrics(c);
    auto t = proxy::time(x1) + c.delay();
    return {proxy::get(x1) != x2, t};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator!= (const T1& x1, const T2& x2)
  {
    return x2 != x1;
  }

  // GREATER THAN
  template<valtype T1, valtype T2>
  val<1> operator> (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) > proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator> (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) > x2, proxy::time(x1)};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator> (const T1& x1, const T2& x2)
  {
    return x2 < x1;
  }
  
  // LESS THAN
  template<valtype T1, valtype T2>
  val<1> operator< (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) < proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator< (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) < x2, proxy::time(x1)};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator< (const T1& x1, const T2& x2)
  {
    return x2 > x1;
  }

  // GREATER THAN OR EQUAL
  template<valtype T1, valtype T2>
  val<1> operator>= (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) >= proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator>= (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) >= x2, proxy::time(x1)};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator>= (const T1& x1, const T2& x2)
  {
    return x2 <= x1;
  }
  
  // LESS THAN OR EQUAL
  template<valtype T1, valtype T2>
  val<1> operator<= (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) <= proxy::get(x2), t};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator<= (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) <= x2, proxy::time(x1)};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator<= (const T1& x1, const T2& x2)
  {
    return x2 >= x1;
  }
  
  // ADDITION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1,T2> operator+ (const T1& x1, const T2& x2)
  {
    static constexpr circuit c = ADD<valt<T1,T2>::size>;
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) + proxy::get(x2), t};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  valt<T1,T2> operator+ (const T1& x1, const T2& x2)
  {
    if (x2==0) return x1;
    static constexpr circuit c = INC<valt<T1>::size>;
    proxy::update_metrics(c);
    auto t = proxy::time(x1) + c.delay();
    return {proxy::get(x1) + proxy::get(x2), t};
  }  

  template<std::integral T1, valtype T2> requires (ival<T2>) // 1st arg constant 
  valt<T1,T2> operator+ (const T1& x1, const T2& x2)
  {
    return x2 + x1;
  }

  // CHANGE SIGN
  template<valtype T>
  valt<T> operator- (const T& x)
  {
    static constexpr circuit c = INV * valt<T>::size + INC<valt<T>::size>;
    proxy::update_metrics(c);
    auto t = proxy::time(x) + c.delay();
    return {-proxy::get(x), t};
  }
  
  // SUBTRACTION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1,T2> operator- (const T1& x1, const T2& x2)
  {
    static constexpr circuit c = SUB<valt<T1,T2>::size>;
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) - proxy::get(x2), t};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  valt<T1,T2> operator- (const T1& x1, const T2& x2)
  {
    if (x2==0) return x1;
    static constexpr circuit c = INC<valt<T1>::size>;
    proxy::update_metrics(c);
    auto t = proxy::time(x1) + c.delay();
    return {proxy::get(x1) - proxy::get(x2), t};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // 1st arg constant 
  valt<T1,T2> operator- (const T1& x1, const T2& x2)
  {
    if (x1==0) return -x2;
    static constexpr circuit c = INC<valt<T2>::size>;
    proxy::update_metrics(c);
    auto t = proxy::time(x2) + c.delay();
    return {proxy::get(x1) - proxy::get(x2), t};
  }
  
  // LEFT SHIFT
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1> operator<< (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) << proxy::get(x2), t};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>)
  valt<T1> operator<< (const T1& x1, const T2& x2)
  {
    // no transistors
    return {proxy::get(x1) << x2, proxy::time(x1)};
  }
  
  // RIGHT SHIFT
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1> operator>> (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) >> proxy::get(x2), t};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>)
  valt<T1> operator>> (const T1& x1, const T2& x2)
  {
    // no transistors
    return {proxy::get(x1) >> proxy::get(x2), proxy::time(x1)};
  }

  // MULTIPLICATION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1,T2> operator* (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) * proxy::get(x2), t};
  }

  template<valtype T1, ival T2> requires (ival<T1>) // second argument is a constant
  valt<T1,T2> operator* (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) * x2, proxy::time(x1)};
  }

  template<ival T1, valtype T2> requires (ival<T2>) // first argument is a constant
  valt<T1,T2> operator* (const T1& x1, const T2& x2)
  {
    return x2 * x1;
  }

  // DIVISION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1> operator/ (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) / proxy::get(x2), t};
  }

  template<valtype T1, ival T2> requires (ival<T1>) // divisor is a constant
  valt<T1> operator/ (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) / x2, proxy::time(x1)};
  }
  
  // MODULO (REMAINDER)
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T2> operator% (const T1& x1, const T2& x2)
  {
    // TODO
    auto t = max(proxy::time(x1),proxy::time(x2));
    return {proxy::get(x1) % proxy::get(x2), t};
  }

  template<valtype T1, ival T2> requires (ival<T1>) // modulus is a constant
  valt<T2> operator% (const T1& x1, const T2& x2)
  {
    // TODO
    return {proxy::get(x1) % x2, proxy::time(x1)};
  }
  
  // BITWISE AND
  template<valtype T1, valtype T2>
  valt<T1,T2> operator& (const T1& x1, const T2& x2)
  {
    static constexpr circuit c = AND<2> * min(T1::size,T2::size);
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) & proxy::get(x2), t};
  }

  template<valtype T1,  std::integral T2> // second argument is constant
  valt<T1> operator& (const T1& x1, const T2& x2)
  {
    // no transistors
    return {proxy::get(x1) & proxy::get(x2), proxy::time(x1)};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator& (const T1& x1, const T2& x2)
  {
    // no transistors
    return x2 & x1;
  }
  
  // BITWISE OR
  template<valtype T1, valtype T2>
  valt<T1,T2> operator| (const T1& x1, const T2& x2)
  {
    static constexpr circuit c = OR<2> * min(T1::size,T2::size);
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) | proxy::get(x2), t};
  }

  template<valtype T1,  std::integral T2> // second argument is constant
  valt<T1> operator| (const T1& x1, const T2& x2)
  {
    // no transistors
    return {proxy::get(x1) | proxy::get(x2), proxy::time(x1)};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator| (const T1& x1, const T2& x2)
  {
    // no transistors
    return x2 | x1;
  }

  // BITWISE EXCLUSIVE OR
  template<valtype T1, valtype T2>
  valt<T1,T2> operator^ (const T1& x1, const T2& x2)
  {
    static constexpr circuit c = XOR<2> * min(T1::size,T2::size);
    proxy::update_metrics(c);
    auto t = max(proxy::time(x1),proxy::time(x2)) + c.delay();
    return {proxy::get(x1) ^ proxy::get(x2), t};
  }

  template<valtype T1,  std::integral T2> // second argument is constant
  valt<T1> operator^ (const T1& x1, const T2& x2)
  {
    circuit c = INV * ones<T1::size>(x2); // not constexpr
    proxy::update_metrics(c);
    return {proxy::get(x1) ^ proxy::get(T1(x2)), proxy::time(x1)+c.delay()};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator^ (const T1& x1, const T2& x2)
  {
    return x2 ^ x1;
  }

  // BITWISE COMPLEMENT
  template<valtype T>
  valt<T> operator~ (const T& x)
  {
    static constexpr circuit c = INV * valt<T>::size;
    proxy::update_metrics(c);
    return {~proxy::get(x), proxy::time(x)+c.delay()};
  }

  // CONCATENATE BITS
  template<valtype T1, valtype T2, valtype... T>
  auto concat(const T1& x1, const T2& x2, const T&... x)
  {
    // no transistors (TODO: wires)
    auto t = std::max({proxy::time(x1),proxy::time(x2),proxy::time(x)...});
    if constexpr (sizeof...(x)==0) {
      auto y = proxy::concatenate(x1,x2);
      proxy::set_time(y,t);
      return y;
    } else {
      auto y = proxy::concatenate(proxy::concatenate(x1,x2),x...);
      proxy::set_time(y,t);
      return y;
    }
  }

  // MAKE ARRAY
  template<typename... T> 
  auto make_array(const T&... args)
  {
    using U = valt<T...>;
    static constexpr u64 N = sizeof...(T);
    arr<U,N> a = {U{args}...};
    return a;
  }
  
  // SELECT BETWEEN TWO VALUES
  template<valtype T>
  valt<T> select(val<1> cond, const T& x1, const T& x2)
  {
    static constexpr auto c = MUX<2,T::size>;
    proxy::update_metrics(c[0]); // MUX select
    proxy::update_metrics(c[1]); // MUX data    
    auto t = std::max({proxy::time(cond)+std::lround(c[0].d),proxy::time(x1),proxy::time(x2)}) + std::lround(c[1].d);
    if (proxy::get(cond)) {
      return {proxy::get(x1), t};
    } else {
      return {proxy::get(x2), t};
    }
  }

  // CONDITIONAL EXECUTION
  // this primitive is too useful not to be provided,
  // however its own delay and energy cost is not modeled, and I do not see an easy way to model it
  template<u64 N, action T>
  void execute(val<N> mask, T f)
  {
    auto prev_exec = exec;
    u64 m = proxy::get(mask);
    auto t = proxy::time(mask);
    for (u64 i=0; i<N; i++) {
      bool cond = (m>>i) & 1;
      exec.set_state(cond,t);
      //if (!cond) continue;
      // to prevent cheating, we execute the action even when the condition is false
      // (otherwise, this primitive could be used to leak any bit)
      if constexpr (std::invocable<T>) {
	f();
      } else {
	f(i);
      }
    }
    exec = prev_exec;
  }
  
  // ###########################
}

#endif // HARCOM_H

