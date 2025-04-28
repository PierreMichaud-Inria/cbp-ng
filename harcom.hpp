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

#ifdef ENABLE_FP
  template<typename T>
  concept arith = std::integral<T> || std::floating_point<T>;
#else
  template<typename T>
  concept arith = std::integral<T>;
#endif

  template<arith auto N>
  struct hard {};

  template<typename T, typename X, typename Y>
  concept unaryfunc = requires (T f, X i) {{f(i)} -> std::convertible_to<Y>;};

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

  template <typename T>
  concept standard_integral = std::integral<T> && ! std::same_as<T,char> && ! std::same_as<T,char8_t> && ! std::same_as<T,char16_t> && ! std::same_as<T,char32_t> && ! std::same_as<T,wchar_t> && ! std::same_as<T,bool>;
  
  template<typename T, typename U>
  bool is_equal(T x, U y)
  {
    if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_equal(x,y);
    } else {
      return x == y;
    }
  }

  template<typename T, typename U>
  bool is_different(T x, U y)
  {
     if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_not_equal(x,y);
    } else {
      return x != y;
    }   
  }
  
  template<typename T, typename U>
  bool is_less(T x, U y)
  {
    if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_less(x,y);
    } else {
      return x < y;
    }
  }

  template<typename T, typename U>
  bool is_greater(T x, U y)
  {
    if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_greater(x,y);
    } else {
      return x > y;
    }    
  }

  template<typename T, typename U>
  bool is_less_equal(T x, U y)
  {
    if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_less_equal(x,y);
    } else {
      return x <= y;
    }
  }

  template<typename T, typename U>
  bool is_greater_equal(T x, U y)
  {
    if constexpr (standard_integral<T> && standard_integral<U>) {
      return std::cmp_greater_equal(x,y);
    } else {
      return x >= y;
    }
  }  

  constexpr auto to_unsigned(std::integral auto x) {return std::make_unsigned_t<decltype(x)>(x);}
  constexpr auto to_unsigned(f32 x) {return std::bit_cast<u32>(x);}
  constexpr auto to_unsigned(f64 x) {return std::bit_cast<u64>(x);}

  template<arith T>
  constexpr u64 bitwidth = (std::same_as<T,bool>)? 1 : sizeof(T)*8;
  
  template<u64 N, std::integral T>
  constexpr auto truncate(T x)
  {
    if constexpr (N >= bitwidth<T>) {
      return x;
    } else {
      return x & ((std::make_unsigned_t<T>(1)<<N)-1);
    }
  }

  template<u64 N, arith T>
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
      y = (y<<1) | (x&1);
    }
    return y;
  }

  constexpr auto reversed_byte = [] () {
    std::array<u8,256> b;
    for (u64 i=0; i<256; i++) {
      b[i] = bit_reversal(i);
    }
    return b;
  }();

  auto constexpr reverse_bits(std::unsigned_integral auto x)
  {
    using T = decltype(x);
    if constexpr (std::same_as<T,bool>) {
      return x;
    } else {
      T y = reversed_byte[x & 0xFF];
      for (u64 i=1; i<sizeof(T); i++) {
	x >>= 8;
	y = (y<<8) | reversed_byte[x & 0xFF];
      }
      return y;
    }
  }

  template<u64 W, std::unsigned_integral T, u64 N>
  auto pack_bits(const std::array<T,N> &in)
  {
    // input array has W-bit elements
    // output array has 64-bit elements
    static_assert(N!=0);
    static constexpr u64 M = (N*W+63)/64;
    if constexpr (W==64) {
      static_assert(std::same_as<T,u64>);
      return in;
    } else {
      static_assert(W!=0 && W<64);
      constexpr u64 mask = (u64(1)<<W)-1;
      std::array<u64,M> out{};
      u64 nbits = 0;
      for (u64 i=0; i<N; i++) {
	u64 x = u64(in[i]) & mask;
	u64 j = nbits / 64;
	u64 pos = nbits % 64;
	assert(j<M);
	out[j] |= x << pos;
	if (pos+W > 64) {
	  assert(pos!=0);
	  assert(j+1<M);
	  out[j+1] = x >> (64-pos);
	}
	nbits += W;
      }
      return out;
    }
  }

  template<u64 W, u64 N>
  auto unpack_bits(const std::array<u64,N> &in)
  {
    // input array has 64-bit elements
    // output array has W-bit elements
    static_assert(N!=0);
    constexpr u64 M = (N*64+W-1)/W;
    if constexpr (W==64) {
      return in;
    } else {
      static_assert(W!=0 && W<64);
      constexpr u64 mask = (u64(1)<<W)-1;
      std::array<u64,M> out{};
      u64 nbits = 0;
      for (u64 j=0; j<M; j++) {
	u64 i = nbits / 64;
	u64 pos = nbits % 64;
	assert(i<N);
	if (pos+W <= 64) {
	  out[j] = (in[i] >> pos) & mask;
	} else {
	  u64 right = in[i] >> pos;	
	  u64 left = (i+1<N)? in[i+1] : 0;
	  out[j] = (right | (left << (64-pos))) & mask;
	}
	nbits += W;
      }
      return out;
    }
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
  // The size of a transistor is limited in practice.
  // However, for consistency with certain limitations and freedoms of the library
  // (no floorplanning except for SRAM, so wire lengths are generally unknown and not modeled,
  // yet the user can define unlimited fanout), we do not limit the maximum scale
  inline constexpr u64 DSMAX = std::numeric_limits<u64>::max(); // default maximum scale is unlimited


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

    constexpr u64 delay() const {return llround(d);}

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
      para.d = std::max(d,x.d);
      para.ci = std::max(ci,x.ci);
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
      para.d = std::max(d,x.d);
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


  // We assume static CMOS implementation for all gates.
  // We do not consider pass-transistor logic (Zimmermann & Fichtner, IEEE JSSC, july 1997)
  // Standard cell libraries use PTL only in a few gates (XOR, MUX)
  // Chinazzo et al., ICECS 2022, "Investigation of pass transistor logic in a 12nm FinFET CMOS technology"

  struct inv : basic_gate { // inverter
    constexpr inv() : basic_gate(2,{INVCAP},INVCAP,INVCAP) {}
  };

  struct nand : basic_gate { // single-input NAND is inverter
    constexpr nand(u64 n) : basic_gate(2*n,{n+GAMMA},n+n*GAMMA,n*(n+GAMMA)) {assert(n>=1);}
  };

  struct nor : basic_gate { // single-input NOR is inverter
    constexpr nor(u64 n) : basic_gate(2*n,{1+n*GAMMA},n+n*GAMMA,n*(1+n*GAMMA)) {assert(n>=1);}
  };

  struct and_nor /*aka AOI21*/: basic_gate { // ~(a+bc)
    constexpr and_nor() : basic_gate(6,{1+2*GAMMA/*a*/,2+2*GAMMA/*b*/,2+2*GAMMA/*c*/},3+2*GAMMA,5+6*GAMMA) {}
  };

  struct or_nand /*aka OAI21*/: basic_gate { // ~(a(b+c))
    constexpr or_nand() : basic_gate(6,{2+GAMMA/*a*/,2+2*GAMMA/*b*/,2+2*GAMMA/*c*/},2+3*GAMMA,6+5*GAMMA) {}
  };

  struct and2_nor /*aka AOI22*/ : basic_gate { // ~(ab+cd)
    constexpr and2_nor() : basic_gate(8,{INVCAP},4+4*GAMMA,8+8*GAMMA) {}
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

  using xnor_cpl = xor_cpl; // a,~a,b,~b ==> ~(a^b) = (~a)^b

  template<f64 STAGE_EFFORT=DSE>
  constexpr u64 num_stages(f64 path_effort, bool odd=false)
  {
    // higher stage effort increases delay but reduces energy
    assert(path_effort>=1);
    i64 n = log(path_effort) / log(STAGE_EFFORT);
    assert(n>=0);
    n += (n & 1) ^ odd;
    return n;
  }

  template<f64 SE=DSE, u64 SMAX=DSMAX>
  constexpr circuit buffer(f64 co, bool cpl, f64 scale=1, f64 bias=0.5)
  {
    f64 ci = inv{}.icap(scale);
    f64 fo = std::min(std::max(1.,co/ci),f64(SMAX));
    u64 ninv = num_stages<SE>(fo,cpl);
    circuit buf;
    if (ninv!=0) {
      f64 foi = pow(fo,1./ninv);
      assert(foi>=1);
      for (u64 i=0; i<ninv; i++) {
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
    // The first segment is driven by a buffer, other segments by large inverters
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
    f64 invscale = std::min(f64(SMAX),optinvscale); // inverter scale
    u64 nseg = std::max(u64(1), u64(llround(length/seglen))); // number of segments
    seglen = length / nseg; // segment length (um)
    f64 cseg = linearcap * seglen; // wire segment capacitance relative to CGATE
    // transistors:
    f64 cfirst = cseg + ((nseg==1)? co : 0);
    circuit c = buffer<SE,SMAX>(cfirst,(nseg-1+cpl)&1,1,bias);
    circuit rep = inv{}.make(cseg,invscale,bias); // repeater = inverter
    for (u64 i=1; i<(nseg-1); i++) c = c + rep;
    if (nseg>1) c = c + inv{}.make(cseg+co,invscale,bias); // last segment drives output cap
    // wire:
    c.w = length;
    c.d += wire_res_delay(METALRES[L]*seglen, cseg*CGATE_pF, rep.ci*CGATE_pF) * (nseg-1);
    c.d += wire_res_delay(METALRES[L]*seglen, cseg*CGATE_pF, co*CGATE_pF); // last segment
    c.e += energy_fJ(METALCAP_fF*length,VDD) * proba_switch(bias);
    if (c.nogate()) c.ci = cfirst; // just a wire
    return c;
  }


  template<u64 B=4/*branching*/>
  constexpr circuit broadcast_tree(u64 n, f64 co=INVCAP, bool cpl=false)
  {
    // broadcast one bit to n outputs with a tree of inverters
    // wires not modeled (TODO)
    static_assert(B>=2);
    if (n<=B) {
      if (cpl) {
	return inv{}.make(co*n);
      } else {
	circuit nogate;
	nogate.ci = co*n;
	return nogate;
      }
    }
    u64 m = n / B;
    u64 r = n % B;
    u64 mm = (r==0)? m : m+1;
    circuit stage = inv{}.make(co*B) * m;
    if (r!=0) stage = stage || inv{}.make(co*r);
    return broadcast_tree<B>(mm,stage.ci,cpl^1) + stage;
  }


  constexpr circuit majority(f64 co)
  {
    // a,b,c ==> ab+ac+bc = ~(~(b+c)+(~a~(bc)))
    and_nor aoi;
    f64 scale_aoi = std::max(1.,sqrt(co/aoi.icap<1>()));
    f64 c1 = aoi.icap<0>(scale_aoi);
    f64 c2 = aoi.icap<1>(scale_aoi);
    circuit i = inv{}.make(c2);
    circuit na = nand{2}.make(c2);
    circuit no = nor{2}.make(c1);
    return ((na | no) || i) + aoi.make(co,scale_aoi); // FIXME: and_nor inputs bias
  }


  template<f64 SE=DSE, u64 SMAX=DSMAX, u64 ARITY=4>
  constexpr circuit nand_nor_tree(u64 n, bool nandfirst, bool cpl, f64 co, f64 scale=1, f64 bias=0.5)
  {
    // alternate NANDs and NORs
    static_assert(ARITY >= 2);
    if (n==0) return {};
    if (n==1) return (cpl)? inv{}.make(co,scale,bias) : circuit{};
    basic_gate gate[2][ARITY+1];
    for (u64 i=1; i<=ARITY; i++) {
      gate[0][i] = nor(i);
      gate[1][i] = nand(i);
    }
    u64 ngates[2][ARITY+1];

    auto partition = [] (u64 n) {
      u64 w = std::min(n,ARITY);
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
    f64 fanout = std::min(co/ci,f64(SMAX));
    u64 path_effort = std::max(1.,path_logical_effort(n)*fanout);
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
      f64 next_scale = std::max(1.,scale * stage_effort / next_le);
      f64 cload = (width!=0)? gate[nand_stage^1][width].icap(next_scale) : (depth>1)? inv{}.icap(next_scale) : co;
      circuit stage;
      for (u64 w=1; w<=ARITY; w++) {
	stage = stage || (gate[nand_stage][w].make(cload,scale,bias) * ngates[nand_stage][w]);
      }
      tree = tree + stage;
      bias = pow(bias,prev_width); // each NAND/NOR stage reduces switching probability (neglect glitching)
      depth--;
      scale = next_scale;
      nand_stage ^= 1;
      if (width!=0) n = std::accumulate(ngates[nand_stage]+1,ngates[nand_stage]+1+ARITY,0);
    } while (width!=0);
    
    if (ninv!=0) {
      // add inverters
      for (u64 i=0; i<ninv; i++) {
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


  constexpr circuit xor2(f64 co, f64 bias=0.5)
  {
    f64 scale = std::max(1.,sqrt(co/(xor_cpl{}.icap())));
    circuit x = xor_cpl{}.make(co,scale,bias);
    circuit i = inv{}.make(x.ci,1.,bias);
    circuit c = i*2+x;
    c.ci = i.ci + x.ci;
    return c;
  }


  constexpr circuit xnor2(f64 co, f64 bias=0.5)
  {
    return xor2(co,bias);
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
	circuit sel = buffer<SE,SMAX>(selcap,false) | buffer<SE,SMAX>(cselcap,true);
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


  template<bool INCR=false>
  constexpr circuit half_adder(f64 co)
  {
    if constexpr (INCR) {
      // a+1
      circuit c = inv{}.make(co);
      c.ci += co;
      return c;
    } else {
       // a+b
      f64 scale = std::max(1.,sqrt(co/(xor_cpl{}.icap()+nor{2}.icap())));
      circuit x = xor_cpl{}.make(co,scale); // sum (a^b)
      circuit n = nor{2}.make(co,scale); // carry out (ab)
      circuit c = x|n;
      circuit i = inv{}.make(c.ci);
      c = i*2 + c;
      c.ci += x.ci;
      return c;
    }
  }


  template<bool INCR=false>
  constexpr circuit full_adder(f64 co)
  {
    if constexpr (INCR) {
      // a+b+1
      f64 scale = std::max(1.,sqrt(co/(xnor_cpl{}.icap()+nand{2}.icap())));
      circuit x = xnor_cpl{}.make(co,scale); // sum = a^~b = ~(a^b)
      circuit n = nand{2}.make(co,scale); // carry = a+b = ~(~a~b)
      circuit xn = x|n;
      circuit c = inv{}.make(xn.ci) * 2 + xn;
      c.ci += xnor_cpl{}.icap(scale);
      return c;
    } else {
      // a+b+c
      // x = ~(abc), y=~(a+b+c), z=ab+ac+bc
      // sum = a^b^c = ~(x(y+z))
      // carry = z
      f64 scale_oai = 1;
      f64 scale_aoi = 1;
      f64 a = and_nor{}.icap<1>() * or_nand{}.icap<1>() / (co*co);
      f64 b = co / or_nand{}.icap<1>();
      auto f = [&](f64 x) {return a*pow(x,4)-b;};
      f64 x0 = sqrt(sqrt(b/a));
      f64 x1 = x0 / (1.-1./(4*a*pow(x0,3)));
      if (x1 >= x0) {
	// always true if co has reasonable value (>=INVCAP)
	assert(f(x0)<=x0 && f(x1)>=x1);
	f64 xm = (x0+x1)/2;
	for (u64 i=0; i<20; i++) {
	  if (f(xm)<xm) {
	    x0 = xm;
	  } else {
	    x1 = xm;
	  }
	}
	scale_aoi = std::max(1.,xm*xm * or_nand{}.icap<1>() / co);
	scale_oai = std::max(1.,xm);
      }
      circuit oai = or_nand{}.make(co,scale_oai);
      f64 coai1 = or_nand{}.icap<0>(scale_oai);
      f64 coai2 = or_nand{}.icap<1>(scale_oai);
      circuit x = nand{3}.make(coai1);
      circuit y = nor{3}.make(coai2);
      circuit aoi = and_nor{}.make(co+coai2,scale_aoi);
      f64 caoi1 = and_nor{}.icap<0>(scale_aoi);
      f64 caoi2 = and_nor{}.icap<1>(scale_aoi);
      circuit na = nand{2}.make(caoi2);
      circuit no = nor{2}.make(caoi1);
      circuit inva = inv{}.make(caoi2);
      circuit z = ((na | no) || inva) + aoi;
      circuit fa = (x|y|z) + oai;
      return fa;
    }
  }


  template<bool INCR=false, bool CARRYIN=false>
  constexpr circuit adder_ks(u64 n, f64 co)
  {
    // Kogge-Stone adder (radix-2)
    // wire capacitance not modeled (TODO?)
    assert(n!=0);
    if (n==1) {
      return (CARRYIN)? full_adder<INCR>(co) : half_adder<INCR>(co);
    }
    u64 depth = std::bit_width(n-1);
    basic_gate G[2] = {or_nand{},and_nor{}}; // inverting generate gate
    basic_gate P[2] = {nor{2},nand{2}}; // inverting propagate gate
    circuit sum = xor2(co) * (n-(CARRYIN^1)); // final sum
    circuit bws = (INCR)? circuit{} : xor2(sum.ci) * n; // bitwise sum
    circuit bwg; // bitwise generate
    circuit bwp; // bitwise propagate
    if constexpr (! INCR) {
      if (n==2) {
	bwg = nand{2}.make(G[0].icap<1>()+sum.ci) * n;
	bwp = nor{2}.make(G[0].icap<1>()+sum.ci) * (n-1);
      } else {
	bwg = nand{2}.make(G[0].icap<0>()+G[0].icap<1>()) * n;
	bwp = nor{2}.make(G[0].icap<1>()+2*P[0].icap()) * (n-1);
      }
    }
    circuit carryout;
    if (depth & 1) {
      carryout.ci = co;
    } else {
      carryout = inv{}.make(co);
    }
    circuit tree = bwp | bwg;
    for (u64 i=0; i<depth; i++) {
      u64 nd = 1<<i; // generate bits already calculated
      u64 ng = n-nd; // G cells
      u64 np = n-std::min(n,2*nd); // P cells
      u64 n2 = n-std::min(n,3*nd); // cells with 2 consumers
      u64 n1 = ng - n2; // cells with single consumer
      u64 np1 = std::min(n1,np); // P cells with single consumer
      if (i < (depth-1)) {
	f64 cog1 = G[(i+1)&1].icap<0>();
	f64 cog2 = G[(i+1)&1].icap<0>() + G[(i+1)&1].icap<1>();
	f64 cop1 = G[(i+1)&1].icap<1>() + P[(i+1)&1].icap();
	f64 cop2 = G[(i+1)&1].icap<1>() + 2 * P[(i+1)&1].icap();
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
    circuit adder = (tree | bws) + (sum || carryout); // large input capacitance (12 CGATE)
    return adder;
  }


  constexpr circuit csa_tree(const std::vector<u64> &icount, f64 co, std::vector<f64> &loadcap)
  {
    // tree of full adders (FA) and half adders (HA)
    // icount vector = number of bits per column (not necessarily uniform)
    // wiring not modeled
    if (icount.empty()) return {};
    u64 n = icount.size(); // number of columns
    assert(n!=0);
    u64 cmax = *std::max_element(icount.begin(),icount.end());
    if (cmax<=2) {
      // final CPA
      assert(loadcap.empty());
      circuit output;
      output.ci = co;
      u64 i2 = 0;
      bool need_cpa = false;
      for (u64 i=0; i<=n; i++) {
	if (i==n || icount.at(i)==0) {
	  if (need_cpa) {
	    circuit cpa = adder_ks(i-i2,co);
	    output = output || cpa;
	    for (u64 j=i2; j<i; j++) loadcap.push_back(cpa.ci);
	    need_cpa = false;
	  }
	  if (i<n) loadcap.push_back(0);
	} else if (! need_cpa) {
	  if (icount.at(i)==1) {
	    loadcap.push_back(co);
	  } else {
	    assert(icount.at(i)==2);
	    need_cpa = true;
	    i2 = i;
	  }
	}
      }
      return output;
    }
    // Wallace tree
    assert(cmax>=3);
    std::vector<u64> n3 (n+1,0);
    std::vector<u64> n2 (n+1,0);
    std::vector<u64> ocount (n+1,0);
    for (u64 i=0; i<n; i++) {
      n3[i] = icount.at(i)/3;
      n2[i] = 0;
      u64 r = icount.at(i)%3;
      if (r==2) {
	n2[i] = 1;
	r = 0;
      }
      ocount.at(i) += n3[i] + n2[i] + r;
      ocount.at(i+1) += n3[i] + n2[i];
    }
    if (ocount.at(n)==0) ocount.pop_back();
    circuit csa = csa_tree(ocount,co,loadcap);
    circuit stage;
    for (u64 i=0; i<n; i++) {
      if (n3[i]==0 && n2[i]==0) continue;
      assert(i+1<loadcap.size());
      f64 ocap = std::max(loadcap[i],loadcap[i+1]);
      circuit col = full_adder(ocap) * n3[i] || half_adder(ocap) * n2[i];
      stage = stage || col;
      loadcap[i] = col.ci;
    }
    return stage + csa;
  }


  constexpr circuit carry_save_adder(const std::vector<u64> &icount, f64 co)
  {
    std::vector<f64> loadcap;
    return csa_tree(icount,co,loadcap);
  }


  constexpr circuit multiplier(u64 n, u64 m, f64 co)
  {
    // multiply two unsigned integers X (n bits) and Y (m bits)
    // X is the multiplier, Y the multiplicand
    if (n>m) {
      return multiplier(m,n,co); // faster and more energy efficient
    }
    assert(m>=n && n>=1);
    u64 cols = m+n-1; // number of CSA columns
    std::vector<u64> count (cols,n); // number of bits per CSA column
    for (u64 i=1; i<n; i++) {
      count[i-1] = i;
      count[cols-i] = i;
    }
    circuit sum = carry_save_adder(count,co); // sum all partial products
    // compute partial products as ~(~a+~b)
    circuit pp1 = nor{2}.make(sum.ci); // partial product (one bit)
    // inputs have high fanout if n or m is large, buffering needed
    circuit bufx = buffer(pp1.ci*m,true); // inverting buffer
    circuit bufy = buffer(pp1.ci*n,true); // inverting buffer
    return (bufx * n || bufy * m) + pp1 * (m*n) + sum;
  }


  template<u64 N>
  constexpr circuit multiplier(u64 m, f64 co)
  {
    // multiply an m-bit unsigned integer multiplicand by a fixed, known multiplier N
    constexpr std::bitset<std::bit_width(N)> n {N};
    u64 cols = m+n.size()-1; // number of CSA columns
    std::vector<u64> count (cols,0); // number of bits per CSA column
    for (u64 i=0; i<n.size(); i++) {
      if (n[i])
	for (u64 j=0; j<m; j++) count.at(i+j)++;
    }
    std::vector<f64> loadcap;
    circuit sum = csa_tree(count,co,loadcap); // sum all partial products
    // inputs may have high fanout, buffering needed
    assert(loadcap.size()>=cols);
    circuit buf;
    for (u64 j=0; j<m; j++) {
      f64 icap = 0;
      for (u64 i=0; i<n.size(); i++) {
	if (n[i]) icap += loadcap.at(i+j);
      }
      buf = buf || buffer(icap,false);
    }
    return buf + sum;
  }


  template<unaryfunc<f64,circuit> F>
  constexpr circuit reduction(u64 n, F op, f64 co)
  {
    // n = number of inputs
    // op represents associative 2-input operation (circuit = op(load_cap))
    // co = output capacitance
    // FIXME: switching activity depends on operation
    assert(n!=0);
    if (n==1) {
      circuit nothing;
      nothing.ci = co;
      return nothing;
    } else if (n==2) {
      return op(co);
    }
    circuit tree = reduction(n-n/2,op,co);
    return op(tree.ci) * (n/2) + tree;
  }


  template<u64 N, unaryfunc<f64,circuit> F>
  constexpr circuit parallel_prefix(F op, f64 (&loadcap)[N])
  {
    // Kogge-Stone parallel prefix tree
    // op represents associative 2-input operation (circuit = op(load_cap))
    // N = number of columns
    // loadcap = output capacitances
    // FIXME: switching activity depends on operation
    static_assert(N!=0);
    circuit tree;
    for (u64 step=std::bit_floor(N-1); step>=1; step/=2) {
      circuit stage;
      for (u64 i=N-1; i>=step; i--) {
	if (i+step*2 < N) loadcap[i] += loadcap[i+step*2];
	circuit c = op(loadcap[i]);
	stage = stage || c;
	loadcap[i] = c.ci;
      }
      tree = stage + tree;
    }
    for (u64 i=0; i<N-1; i++)
      loadcap[i] += loadcap[i+1];
    tree.ci = *std::max_element(loadcap,loadcap+N);
    return tree;
  }


  template<u64 N>
  constexpr circuit parallel_prefix(u64 data, circuit (*op1)(f64), circuit (*op2)(f64), f64 (&loadcap)[N])
  {
    // Kogge-Stone parallel prefix tree but with alternating polarity:
    // 1st stage uses op1, 2nd stage uses op2, 3rd stage uses op1, 4th stage uses op2, and so on
    // The actual associative operation is A op B = ~(A op1 B) = (~A op2 ~B)
    // For example, op1=NOR implies op2=NAND and op=OR
    // data = number of data bits
    // N = number of columns
    // loadcap = output capacitances
    // FIXME: switching activity depends on operation
    static_assert(N!=0);
    std::array op = {op1,op2};
    u64 stages = std::bit_width(N-1); // number of computation stages
    bool pol = stages & 1; // polarity of the output of the last computation stage (pol=1 means inverted)
    circuit tree;
    if (pol) {
      // odd number of computation stages: end with a stage of inverters
      for (u64 i=0; i<N; i++) {
	circuit c = inv{}.make(loadcap[i]) * data;
	tree = tree || c;
	loadcap[i] = c.ci;
      }
    }
    for (u64 step=std::bit_floor(N-1); step>=1; step/=2) {
      circuit stage;
      for (u64 i=N-1; i>=step; i--) {
	if (i+step*2 < N) loadcap[i] += loadcap[i+step*2];
	circuit c = op[pol](loadcap[i]);
	stage = stage || c;
	loadcap[i] = c.ci;
      }
      // inverters are used to propagate values that are computed earlier than the last stage
      // (FIXME: some inverters are unnecessary)
      for (u64 i=0; i<step; i++) {
	circuit c = inv{}.make(loadcap[i]) * data;
	stage = stage || c;
	loadcap[i] = c.ci;
      }
      tree = stage + tree;
      pol ^= 1;
    }
    assert(!pol);
    for (u64 i=0; i<N-1; i++)
      loadcap[i] += loadcap[i+1];
    tree.ci = *std::max_element(loadcap,loadcap+N);
    return tree;    
  }


  constexpr circuit unsigned_greater_than(u64 n, f64 co, bool gte = false)
  {
    // Inputs: two n-bit unsigned integers A=(An,...,A1) and B=(Bn,...,B1)
    // Outputs: G=(A>B), E=(A==B)
    // Divide and conquer: A=XY, B=ZW ==> G = (X>Z) || (X==Z && Y>W)
    // This 2-bit binary operation is associative: (g,e) * (g',e') = (g+eg',ee')
    // Define (G[i:i],E[i:i]) = (Ai~Bi,~(Ai^Bi)) for all i in [1,n]
    // Define, for i>j>k, (G[i:k],E[i:k]) = (G[i:j],E[i:j]) * (G[j-1:k],E[j-1:k])
    // We have: G=G[n:1] and E=E[n:1]
    // Use reduction tree with alternating polarity: (AOI,NAND) at stage 1, (OAI,NOR) at stage 2
    // NB: A>=B is equivalent to ~(A<B) ==> start with inverted polarity
    assert(n!=0);
    basic_gate I[2] = {nor{2},nand{2}}; // initial G stage, depending on polarity
    basic_gate G[2] = {and_nor{},or_nand{}}; // reduction stage G, depending on polarity
    basic_gate E[2] = {nand{2},nor{2}}; // reduction stage E, depending on polarity      
    if (n==1) {
      circuit c = inv{}.make(I[gte].icap()) + I[gte].make(co);
      c.ci = std::max(c.ci,I[gte].icap());
      return c;
    }
    assert(n>=2);
    bool pol = gte; // polarity of the inputs (pol=1 means inverted)
    // initial stage
    f64 cgleft = G[pol].icap<0>();
    f64 celeft = G[pol].icap<1>() + E[pol^1].icap();
    f64 cgright = G[pol].icap<2>();
    f64 ceright = E[pol].icap();          
    circuit gleft = (inv{}.make(I[gte].icap()) + I[gte].make(cgleft)) * (n/2);
    gleft.ci = std::max(gleft.ci,I[gte].icap());
    circuit gright = (inv{}.make(I[gte].icap()) + I[gte].make(cgright)) * (n/2);
    gright.ci = std::max(gright.ci,I[gte].icap());
    circuit ginv = (inv{}.make(I[gte].icap()) + I[gte].make(inv{}.icap())) * (n%2);
    ginv.ci = std::max(ginv.ci,I[gte].icap());
    circuit eleft = xnor2(celeft) * (n/2); // XOR and XNOR are identical circuits
    circuit eright = xnor2(ceright) * (n/2-1); // the rightmost E circuit is not needed
    circuit einv = xnor2(inv{}.icap()) * (n%2);
    circuit tree = (gleft | eleft) || (gright | eright) || (ginv | einv);
    // reduction tree
    while (n>=2) {
      u64 nops = n/2; // number of reduction operations in this stage
      u64 nextn = n-nops; // value of n for the next stage
      circuit stage;
      if (nextn == 1) {
	// last reduction stage
	f64 ocap = (pol==gte)? inv{}.icap() : co;
	stage = G[pol].make(ocap); // the E circuit is not needed
      } else {
	// next stage is a reduction stage
	f64 cgleft = G[pol^1].icap<0>();
	f64 celeft = G[pol^1].icap<1>() + E[pol^1].icap();
	f64 cgright = G[pol^1].icap<2>();
	f64 ceright = E[pol^1].icap();
	u64 nextnops = nextn/2; // >=1
	u64 nopsfeedright = nextnops; // >=1
	u64 nopsfeedleft = nextnops - (n%2) * (1-nextn%2);
	// when n is odd, transmit the remaining bundle (g,e) through a pair of inverters
	u64 nopsfeedinv = (1-n%2) * (nextn%2);
	u64 ninvfeedleft = nextnops - nopsfeedleft;
	u64 ninvfeedinv = (n%2) * (nextn%2);
	circuit opsfeedleft = (G[pol].make(cgleft) | E[pol].make(celeft)) * nopsfeedleft;
	// the rightmost E circuit is not needed
	circuit opsfeedright = G[pol].make(cgright) * nopsfeedright | E[pol].make(ceright) * (nopsfeedright-1);
	circuit opsfeedinv = (G[pol].make(INVCAP) | E[pol].make(INVCAP)) * nopsfeedinv;
	circuit invfeedleft = (inv{}.make(cgleft) | inv{}.make(celeft)) * ninvfeedleft;
	circuit invfeedinv = inv{}.make(INVCAP) * (2 * ninvfeedinv); // FIXME: unneeded
	stage = opsfeedleft || opsfeedright || opsfeedinv || invfeedleft || invfeedinv;
      }
      tree = tree + stage;
      n = nextn;
      pol ^= 1;
    }
    // pol is the polarity of the outputs of the reduction tree
    if (pol != gte) {
      tree = tree + inv{}.make(co);
    }
    return tree;
  }


  // ###########################

  template<u64 N>
  struct flipflops {
    static constexpr circuit oneflop = []() {
      // master-slave flip-flop (2 latches)
      inv i;
      inv_tri t;
      f64 co = INVCAP; // we do not care about the delay anyway
      f64 cfeedback1 = i.icap() + t.cp;
      f64 cfeedback2 = 2*i.icap() + t.cp;
      circuit tri1 = t.make(cfeedback1);
      circuit tri2 = t.make(cfeedback2);
      circuit forward = tri1 + i.make(2*t.icap()) + tri2 + i.make(co);
      circuit backward = tri1 || (i.make(t.icap()) + tri2);
      return forward || backward;
    } ();

    static constexpr circuit flops = oneflop * N;

    // master-latch clock generated locally at each flop
    static constexpr circuit clock1 = inv{}.make(inv_tri{}.icap<1>()*4);

    static constexpr f64 width = SRAM_CELL.wordline_length * N; // um
    static constexpr f64 height = SRAM_CELL.bitline_length * (oneflop.t+clock1.t) / SRAM_CELL.transistors; // um

    static constexpr circuit clocking = []() {
      f64 phicap2 = inv_tri{}.icap<2>() * 4;
      circuit clock2 = wire<4.>(width,false,(phicap2+clock1.ci)*N); // slave-latch clock
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
    static constexpr u64 XTORS = T::num_xtors();
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
    static constexpr u64 num_xtors() {return 0;}
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
    static_assert(M%D==0);
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
    static constexpr f64 SACAP = std::max(SACAPMIN,std::min(SACAPMAX,N*BLCAP*XCSA)); // sense amp capacitance relative to CGATE
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
    static constexpr f64 array_width() {return sqrt(AREA*ASPECT_RATIO);} // um
    static constexpr f64 array_height() {return sqrt(AREA/ASPECT_RATIO);} // um
    static constexpr u64 num_xtors() {return CELL_XTORS + PERI_XTORS;}

    static constexpr f64 read_latency()
    {
      return std::max(WLSEL.d + WLRC + BLDELAY + SADELAY, CMUX[0].d) + CMUX[1].d; // ps
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
    static constexpr f64 read_latency() { return ABUS.d + BANK::LATENCY;}
    static constexpr f64 read_energy() {return ABUS.e + BANK::EREAD * NB + BANKR::EREAD;}
    static constexpr f64 write_energy() {return ABUS.e + BANK::EWRITE * NB + BANKR::EWRITE;}
    static constexpr f64 array_width() {return BANK::WIDTH * NB + BANKR::WIDTH;}
    static constexpr f64 array_height() {return BANK::HEIGHT;}
    static constexpr u64 num_xtors() {return BANK::XTORS * NB + BANKR::XTORS + ABUS.t;}

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
    static constexpr f64 read_latency() {return READ.d;}
    static constexpr f64 read_energy() {return READ.e;}
    static constexpr f64 write_energy() {return ACC.e + WTREE.e + BANK::EWRITE;}
    static constexpr f64 array_width() {return BANK::WIDTH * BX;}
    static constexpr f64 array_height() {return BANK::HEIGHT * BY;}

    static constexpr u64 num_xtors()
    {
      return NB * BANK::XTORS + ACC.t + WTREE.t + RTREE[0].t + RTREE[1].t;
    }
  };


  constexpr bool ok_config(u64 E/*entries*/, u64 D/*data bits*/, u64 MAXN, u64 MAXM)
  {
    bool ok = E!=0 && D!=0; // mandatory
    ok = ok && std::has_single_bit(MAXN); // mandatory: MAXN must be a power of 2
    // prune off configs that are likely bad
    u64 banksize = MAXN * std::max(D,MAXM);
    ok = ok && (16*E*D <= banksize*banksize || (MAXN>=1024 && MAXM>=512));
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
      constexpr f64 FACTOR = std::max(2.,BIAR*BIAR);
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
	constexpr u64 BY = std::min(to_pow2(S/N),std::bit_floor(R/N));
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
    static constexpr u64 num_xtors() {return ARR::XTORS;}
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
  inline constexpr circuit XOR = reduction(N,[](f64 co){return xor2(co);},OUTCAP);

  template<u64 WIDTH>
  inline constexpr circuit ADD = adder_ks<false,false>(WIDTH,OUTCAP);

  template<u64 WIDTH>
  inline constexpr circuit INC = adder_ks<true,true>(WIDTH,OUTCAP);

  template<u64 WIDTH>
  inline constexpr circuit SUB = adder_ks<false,true>(WIDTH,OUTCAP);

  template<u64 WIDTH, u64 N> requires (N>=2)
  inline constexpr circuit ADDN = []() {
    // add N integers, WIDTH bits each
    std::vector<u64> count (WIDTH,N);
    return carry_save_adder(count,OUTCAP);
  }();

  template<u64 N>
  inline constexpr circuit EQUAL = xnor2(AND<N>.ci) * N + AND<N>;

  template<u64 N>
  inline constexpr circuit NEQ = xor2(OR<N>.ci) * N + OR<N>;

  template<u64 N>
  inline constexpr circuit GT = unsigned_greater_than(N,OUTCAP,false); // TODO: signed

  template<u64 N>
  inline constexpr circuit GTE = unsigned_greater_than(N,OUTCAP,true); // TODO: signed

  template<u64 N, u64 WIDTH> requires (N>=2)
  inline constexpr auto MUX = mux(N,WIDTH,OUTCAP);

  template<u64 N>
  constexpr circuit priority_encoder = []() {
    static_assert(N!=0);
    if constexpr (N==1) {
      return circuit{};
    } else {
      f64 co = INVCAP;
      circuit output = nor{2}.make(co);
      f64 ocap[N-1];
      std::fill(ocap,ocap+N-1,output.ci);
      auto nor2 = [] (f64 co) {return nor{2}.make(co);};
      auto nand2 = [] (f64 co) {return nand{2}.make(co);};
      circuit prefix_or = parallel_prefix(1,nor2,nand2,ocap);
      circuit input = inv{}.make(prefix_or.ci);
      return input * (N-1) + prefix_or + output * (N-1);
    }
  }();  

  // circuit for replicating a value M times
  template<u64 DATABITS, u64 COPIES>
  inline constexpr circuit REP = []() {
    // wires not modeled (TODO)
    circuit tree = broadcast_tree(COPIES);
    if (tree.ci > INVCAP) {
      // composing replication must not allow to build a large replicator at no hardware cost,
      // so we charge a hardware cost corresponding to that of an inverter
      tree = inv{}.make(tree.ci) + tree;
    }
    return tree * DATABITS;
  }();

  // signed mul has roughly same complexity as unsigned
  template<u64 N, u64 M>
  inline constexpr circuit IMUL = multiplier(N,M,OUTCAP); // N-bit x M-bit

  template<u64 HARDN, u64 M>
  inline constexpr circuit HIMUL = multiplier<HARDN>(M,OUTCAP); // M-bit x HARDN


  // ###########################

  template<typename T>
  struct toval_impl {};

  template<arith T>
  struct toval_impl<T> {
    using type = val<bitwidth<T>,T>;
  };

  template<valtype T>
  struct toval_impl<T> {
    using U = std::remove_reference_t<T>;
    using type = val<U::size,typename U::type>;
  };

  template<typename T>
  using toval = toval_impl<T>::type;

  template<typename T>
  struct base_impl {};

  template<arith T>
  struct base_impl<T> {using type = T;};
  
  template<valtype T>
  struct base_impl<T> {using type = toval<T>::type;};

  template<typename T>
  using base = base_impl<T>::type;

  template<typename T>
  concept ival = std::integral<std::remove_reference_t<T>> || std::integral<base<T>>;

  template<typename T>
  concept fval = std::floating_point<std::remove_reference<T>> || std::floating_point<base<T>>;

  template<typename T>
  inline constexpr u64 length = 0;

  template<valtype T>
  inline constexpr u64 length<T> = toval<T>::size;

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

  template<action A>
  struct action_return {};
  
  template<action A> requires (std::invocable<A>)
  struct action_return<A> {
    using type = std::invoke_result_t<A>;
  };

  template<action A> requires (std::invocable<A,u64>)
  struct action_return<A> {
    using type = std::invoke_result_t<A,u64>;
  };

  template<action A>
  using return_type = action_return<A>::type;


  // ###########################
  // GLOBAL VARIABLES

  inline class exec_control {
    template<u64,arith> friend class val;
    template<memdatatype,u64> friend class ram;
    friend class globals;
    template<valtype T, action A> requires (std::same_as<return_type<A>,void>)
      friend void execute(T&&,const A&);
    template<valtype T, action A> friend auto execute(T&&,const A&);
  private:
    bool active = true;
    u64 time = 0;
    exec_control(const exec_control &s) = default;
    exec_control& operator=(const exec_control &s) = default;
    void set_state(bool cond, u64 t=0) {active = cond; time = t;}
  public:
    exec_control() : active(true), time(0) {}
  } exec;


  template<typename T> requires (std::integral<T> || std::floating_point<T>)
  class global {
    friend class globals;
    friend class ::simulator;
    template<u64,arith> friend class val;
    template<u64,arith> friend class reg;
    template<memdatatype,u64> friend class ram;
  private:
    T data = 0;
    global() : data(0) {}
    global(const global&) = default;
    global& operator= (const global& x) = default;
    void operator+= (T i) {data+=i;}
    void operator++(int) {data++;}
    void operator&() = delete;
  public:
    global(T x) : data(x) {}
    operator T() const {return data;}

    void print(std::string s = "", std::ostream & os = std::cout) const
    {
      os << s << data << std::endl;
    }
  };


  inline class globals {
    template<u64,arith> friend class val;
    template<u64,arith> friend class reg;
    template<valtype,u64> friend class arr;
    template<memdatatype,u64> friend class ram;
    friend class proxy;
    friend class ::simulator;
    static constexpr u64 first_cycle = 1;
  public:
    global<u64> clock_cycle_ps = 300;
    global<u64> cycle = first_cycle;
    global<u64> storage;
    global<f64> energy_fJ;
    global<u64> transistors_storage;
    global<u64> transistors_logic[2];

  private:
    global<bool> storage_destroyed = false;

    void next_cycle()
    {
      assert(clock_cycle_ps != 0);
      cycle++;
      transistors_logic[1] = transistors_logic[0];
      transistors_logic[0] = 0;
    }

    void update_storage(u64 nbits, u64 xtors)
    {
      storage += nbits;
      transistors_storage += xtors;
    }

    void update_energy(f64 e)
    {
      if (exec.active) energy_fJ += e;
    }

    void update_logic(const circuit &c)
    {
      transistors_logic[0] += c.t;
      update_energy(c.e);
    }

  public:
    f64 power_mW()
    {
      assert(cycle>first_cycle);
      assert(clock_cycle_ps != 0);
      return energy_fJ / ((cycle-first_cycle)*clock_cycle_ps);
    }

    void print(std::ostream & os = std::cout)
    {
      os << std::setprecision(3);
      if (clock_cycle_ps != 0 && cycle>first_cycle) {
	clock_cycle_ps.print("clock cycle (ps): ",os);
	os << "clock frequency (GHz): " << 1000./clock_cycle_ps << std::endl;
      }
      storage.print("storage (bits): ",os);
      transistors_storage.print("transistors (storage): ",os);
      if (cycle == first_cycle) {
	transistors_logic[0].print("transistors (logic): ",os);
      } else {
	transistors_logic[1].print("transistors (logic): ",os);
      }
      if (cycle == first_cycle) {
	energy_fJ.print("energy (fJ): ",os);	
      } else if (clock_cycle_ps != 0) {
	os << "power (mW): " << power_mW() << std::endl;
      }
    }
  } panel;

  
  // ###########################


  template<u64 N, arith T = u64>
  class val {
    //static_assert(N!=0);
    static_assert(N<=bitwidth<T>,"the number of bits exceeds the underlying C++ type");
    static_assert(N==bitwidth<T> || std::integral<T>);

    template<u64,arith> friend class val;
    template<u64,arith> friend class reg;
    template<u64,u64> friend class split;
    template<valtype, u64> friend class arr;
    template<memdatatype,u64> friend class ram;
    template<u64,u64,arith> friend class rom;
    friend class proxy;
    friend class ::simulator;

  private:
    T data = 0;
    u64 timing = 0; // time
    u64 read_credit = 0;

    T fit(T x) const requires std::integral<T>
    {
      return truncate<N>(x);
    }

    T sign_extended() const requires std::signed_integral<T>
    {
      if constexpr (N < bitwidth<T>) {
	assert((std::make_unsigned_t<T>(data) >> N) == 0);
      }
      constexpr T signbit = T(1) << (N-1);
      return (data ^ signbit) - signbit;
    }

    void set_time(u64 t)
    {
      timing = t;
    }

    void read()
    {
      // simulation speed is impacted by this function
#ifndef FREE_FANOUT
      if (read_credit!=0) {
	read_credit--;
      } else {
#ifdef CHECK_FANOUT
	assert(read_credit!=0 && "fanout exhausted");
#else
	// fanout exhausted: delay increases linearly with the number of reads
	// delay increment is that of a FO2 inverter (wires not modeled, TODO?)	
	constexpr circuit fo2inv = inv{}.make(2*INVCAP) * N; // N inverters in parallel
	static_assert(fo2inv.delay()!=0);
	panel.update_logic(fo2inv);
	set_time(time()+fo2inv.delay());
#endif // CHECK_FANOUT
      }
#endif // FREE_FANOUT
    }

    T get() & // lvalue
    {
      read();
      return data;
    }

    T get() && // rvalue
    {
      auto old_data = data;
      data = 0; // destructive read
      return old_data;
    }

    T get() & requires std::signed_integral<T> // lvalue
    {
      read();
      return sign_extended();
    }

    T get() && requires std::signed_integral<T> // rvalue
    {
      auto old_data = sign_extended();
      data = 0; // destructive read
      return old_data;
    }

    auto get_vt() & // lvalue
    {
      return std::tuple {get(),time()}; // list initialization, get() executes before time()
    }

    auto get_vt() && // rvalue
    {
      return std::tuple {std::move(*this).get(),time()}; // list initialization
    }    

    void operator= (val & x)
    {
      auto [xv,xt] = x.get_vt();
      if (exec.active) data = xv;
      set_time(xt);
    }

    void operator= (val && x)
    {
      auto [xv,xt] = std::move(x).get_vt();
      if (exec.active) data = xv;
      set_time(xt);
    }

    void operator&() = delete;

  public:

    static constexpr u64 size = N;
    using type = T;

    u64 time() const
    {
      return std::max(exec.time,timing);
    }

    static constexpr T maxval = []() {
      if constexpr (N == bitwidth<T>) {
	return std::numeric_limits<T>::max();
      } else {
	static_assert(N < bitwidth<T>);
	if constexpr (std::unsigned_integral<T>) {
	  return (T(1)<<N)-1;
	} else {
	  static_assert(std::signed_integral<T>);
	  return (T(1)<<(N-1))-1;
	}
      } 
    }();

    static constexpr T minval = []() {
      if constexpr (N == bitwidth<T>) {
	return std::numeric_limits<T>::min();
      } else {
	static_assert(N < bitwidth<T>);
	if constexpr (std::unsigned_integral<T>) {
	  return 0;
	} else {
	  static_assert(std::signed_integral<T>);
	  return -(T(1)<<(N-1));
	}
      } 
    }();

    val() {}

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
      static_assert(N == bitwidth<T>);
      set_time(t);
    }

    val(val &x) : val{x.get(),x.time()} {} // list initialization, get() executes before time()

    val(const val &&x) : val{std::move(x).get(),x.time()} {} // list initialization

    template<valtype U> requires std::unsigned_integral<T>
    val(U &&x) : val{to_unsigned(std::forward<U>(x).get()),x.time()} {} // list initialization

    template<valtype U>
    val(U &&x) : val{std::forward<U>(x).get(),x.time()} {} // list initialization

    template<arith auto FO>
    void fanout(hard<FO>) & // lvalue
    {
#ifndef FREE_FANOUT
      static_assert(FO>=2);
      // delay logarithmic with fanout (wires not modeled, TODO?)
      panel.update_logic(REP<N,FO>);
      set_time(time()+REP<N,FO>.delay());
      read_credit = FO;
#endif
    }

    [[nodiscard]] val&& fo1() & // lvalue
    {
      return std::move(*this);
    }

    void fanout(auto) && = delete; // rvalue, not needed

    void print(std::string before="", std::string after="\n", bool t=true, std::ostream & os=std::cout) const
    {
      if constexpr (std::signed_integral<T>) {
	os << before << +sign_extended();
      } else {
	os << before << +data;
      }
      if (t)
	os << " (t=" << time() << " ps)";
      os << after << std::flush;
    }

    void printb(std::string before="", std::string after="\n", bool t=true, std::ostream & os=std::cout) const
    {
      os << before;
      if constexpr (std::integral<T>) {
	os << std::bitset<N>(data);
      } else if constexpr (std::same_as<T,f32>) {
	os << std::bitset<N>(std::bit_cast<u32>(data));
      } else if constexpr (std::same_as<T,f64>) {
	os << std::bitset<N>(std::bit_cast<u64>(data));
      }
      if (t)
	os << " (t=" << time() << " ps)";
      os << after << std::flush;
    }

    template<u64 W>
    auto make_array(val<W>&&) & // lvalue
    {
      static_assert(N!=0);
      static_assert(std::unsigned_integral<T>);
      static_assert(W!=0 && W<=64);
      constexpr u64 mask = (W==64)? -1 : (u64(1)<<W)-1;
      constexpr u64 M = (N+W-1)/W;
      auto [x,t] = get_vt();
      arr<val<W>,M> out;
      for (u64 i=0; i<M; i++) {
	out[i] = val<W>{x&mask,t};
	x >>= W;
      }
      return out;
    }

    template<u64 W>
    auto make_array(val<W>&&) && // rvalue
    {
      static_assert(N!=0);
      static_assert(std::unsigned_integral<T>);
      static_assert(W!=0 && W<=64);
      constexpr u64 mask = (W==64)? -1 : (u64(1)<<W)-1;
      constexpr u64 M = (N+W-1)/W;
      auto [x,t] = std::move(*this).get_vt();
      arr<val<W>,M> out;
      for (u64 i=0; i<M; i++) {
	out[i] = val<W>{x&mask,t};
	x >>= W;
      }
      return out;
    }

    [[nodiscard]] val reverse() & requires std::unsigned_integral<T> // lvalue
    {
      // no transistors
      return {reverse_bits(get()) >> (bitwidth<T>-N), time()};
    }

    [[nodiscard]] val reverse() && requires std::unsigned_integral<T> // rvalue
    {
      // no transistors
      return {reverse_bits(std::move(*this).get()) >> (bitwidth<T>-N), time()};
    }

    auto ones() & // lvalue
    {
      constexpr circuit c = (N>=2)? ADDN<1,N> : circuit{};
      panel.update_logic(c);
      auto [x,t] = get_vt();
      auto n = std::popcount(truncate<N>(x));
      return val<std::bit_width(N)> {n,t+c.delay()};
    }

    auto ones() && // rvalue
    {
      constexpr circuit c = (N>=2)? ADDN<1,N> : circuit{};
      panel.update_logic(c);
      auto [x,t] = std::move(*this).get_vt();
      auto n = std::popcount(truncate<N>(x));
      return val<std::bit_width(N)> {n,t+c.delay()};
    }

    [[nodiscard]] val priority_encode() & requires std::unsigned_integral<T> // lvalue
    {
      constexpr circuit c = priority_encoder<N>;
      auto [x,t] = get_vt();
      u64 y = x & (x^(x-1));
      panel.update_logic(c);
      return {y,t+c.delay()};
    }

    [[nodiscard]] val priority_encode() && requires std::unsigned_integral<T> // rvalue
    {
      constexpr circuit c = priority_encoder<N>;
      auto [x,t] = std::move(*this).get_vt();
      u64 y = x & (x^(x-1));
      panel.update_logic(c);
      return {y,t+c.delay()};
    }

    template<arith auto M>
    arr<val,M> replicate(hard<M>) & // lvalue
    {
      // only the user knows the actual fanout (>=M) and can set it
      static_assert(M>=2);
      arr<val,M> a;
      for (u64 i=0; i<M; i++) a[i] = *this;
      return a;
    }

    template<arith auto M>
    arr<val,M> replicate(hard<M>) && // rvalue
    {
      // the user cannot set the fanout (rvalue), but the fanout is known
      static_assert(M>=2);
      fanout(hard<M>{});
      arr<val,M> a;
      for (u64 i=0; i<M; i++) a[i] = std::move(*this);
      return a;
    }

  };


  // ###########################

  class proxy {
  private:
    proxy() = delete;
    ~proxy() = delete;

    template<typename T>
    static auto get(T &&x)
    {
      if constexpr (valtype<T> || arrtype<T>) {
	return std::forward<T>(x).get();
      } else {
	static_assert(arith<std::remove_reference_t<T>>);
	return x;
      }
    }

    template<typename T>
    static u64 time(T &&x)
    {
      if constexpr (valtype<T> || arrtype<T>) {
	return std::forward<T>(x).time();
      } else {
	static_assert(arith<std::remove_reference_t<T>>);
	return 0;
      }
    }

    template<typename T>
    static auto get_vt(T &&x)
    {
      if constexpr (valtype<T>) {
	return std::forward<T>(x).get_vt();
      } else {
	using atype = std::remove_reference_t<T>;
	static_assert(arith<atype>);
	return std::tuple<atype,u64>{x,0};
      }
    }

    template<typename T>
    static void set_time(T &&x, u64 t)
    {
      static_assert(valtype<T> || arrtype<T>);
      std::forward<T>(x).set_time(t);
    }

    static void update_logic(const circuit &c)
    {
      panel.update_logic(c);
    }

    template<valtype T1, valtype T2>
    static auto concat2(T1 &&x1, T2 &&x2)
    {
      static_assert(std::unsigned_integral<base<T1>>);
      static_assert(std::unsigned_integral<base<T2>>);
      using rtype = val<x1.size+x2.size>;
      if constexpr (x1.size == 0) {
	auto [v2,t2] = std::forward<T2>(x2).get_vt();
	return rtype {v2,t2};
      } else {
	auto [v1,t1] = std::forward<T1>(x1).get_vt();
	auto [v2,t2] = std::forward<T2>(x2).get_vt();
	return rtype {(v1 << x2.size) | v2, std::max(t1,t2)};
      }
    }

    template<u64,arith> friend class reg;
    template<valtype,u64> friend class arr;
    template<u64,u64> friend class split;
    template<memdatatype,u64> friend class ram;
    template<u64,u64,arith> friend class rom;
    
    template<valtype T1, valtype T2>
    friend val<1> operator== (T1&&, T2&&);

    template<valtype T1, arith T2>
    friend val<1> operator== (T1&&, T2);
    
    template<valtype T1, valtype T2>
    friend val<1> operator!= (T1&&, T2&&);

    template<valtype T1, arith T2>
    friend val<1> operator!= (T1&&, T2);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend val<1> operator> (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend val<1> operator> (T1&&, T2);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend val<1> operator< (T1&&, T2&&);
    
    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend val<1> operator< (T1&&, T2);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend val<1> operator>= (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend val<1> operator>= (T1&&, T2);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend val<1> operator<= (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend val<1> operator<= (T1&&, T2);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend auto operator+ (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend auto operator+ (T1&&, T2);

    template<valtype T>
    friend auto operator- (T&&);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend auto operator- (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend auto operator- (T1&&, T2);
    
    template<std::integral T1, valtype T2> requires (ival<T2>)
    friend auto operator- (T1, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1> operator<< (T1&&, T2);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1> operator>> (T1&&, T2);
    
    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend auto operator* (T1&&, T2&&);

    template<valtype T1, arith auto N2> requires (ival<T1>)
    friend auto operator* (T1&&, hard<N2>);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T1> operator/ (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T1> operator/ (T1&&, T2);

    template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
    friend valt<T2> operator% (T1&&, T2&&);

    template<valtype T1, std::integral T2> requires (ival<T1>)
    friend valt<T2> operator% (T1&&, T2);
    
    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator& (T1&&, T2&&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator& (T1&&, T2);

    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator| (T1&&, T2&&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator| (T1&&, T2);
    
    template<valtype T1, valtype T2>
    friend valt<T1,T2> operator^ (T1&&, T2&&);

    template<valtype T1, std::integral T2>
    friend valt<T1> operator^ (T1&&, T2);
    
    template<valtype T>
    friend valt<T> operator~ (T&&);

    template<valtype T1, valtype T2, valtype... T>
    friend auto concat(T1&&, T2&&, T&&...);

    template<valtype T, valtype T1, valtype T2>
    friend valt<T1,T2> select(T&&, T1&&, T2&&);

    template<valtype T, action A> requires (std::same_as<return_type<A>,void>)
    friend void execute(T&&,const A&);

    template<valtype T, action A>
    friend auto execute(T&&,const A&);  
  };


  // ###########################

  template<u64 L, u64 R>
  class split {
  public:
    val<L> left;
    val<R> right;

    template<valtype T>
    split(T && x)
    {
      static_assert(x.size == L+R);
      auto [data,t] = std::forward<T>(x).get_vt();
      right = {data,t};
      if constexpr (R<64) {
	left = {data>>R,t};
      }
    }
  };


  // ###########################


  template<u64 N, arith T = u64>
  class reg : public val<N,T> {
  public:
    using stg = flipflops<N>;

  private:
    u64 last_write_cycle = 0;

    void create()
    {
      assert(!panel.storage_destroyed && "all storage (reg,ram) must have the same lifetime");
      panel.update_storage(N,stg::xtors);
    }

  public:

    reg()
    {
      create();
    }

    reg(reg &other) : val<N,T>{other}
    {
      create();
    }

    reg(reg &&) = delete;

    template<std::convertible_to<val<N,T>> U>
    reg(U &&x) : val<N,T>{std::forward<U>(x)}
    {
      create();
    }    

    ~reg()
    {
      panel.storage_destroyed = true;
    }

    template<typename U>
    void assign(U &&other)
    {
      assert(panel.cycle>last_write_cycle && "single register write per cycle");
      last_write_cycle = panel.cycle;     
      val<N,T>::operator=(std::forward<U>(other));
      panel.update_energy(stg::write_energy_fJ);
    }

    void operator= (reg & other)
    {
      assign(other);
    }

    void operator= (reg &&) = delete;

    template<typename U>
    void operator= (U &&other)
    {
      assign(std::forward<U>(other));
    }

  };


  // ###########################

  template<valtype T, u64 N>
  class arr {
    template<valtype,u64> friend class arr;
    template<u64,arith> friend class val;
    template<memdatatype,u64> friend class ram;
    friend class proxy;
  public:
    static constexpr u64 size = N;
    using type = T;
    using atype = T::type;

  private:

    std::array<T,N> elem {};

    template<arrtype U>
    void copy_from(U &&other)
    {
      static_assert(other.size == N);
      if constexpr (std::is_rvalue_reference_v<decltype(other)>) {
	for (u64 i=0; i<N; i++) elem[i] = std::move(other.elem[i]);
      } else {
	for (u64 i=0; i<N; i++) elem[i] = other.elem[i];
      }
    }

    void operator= (arr &other)
    {
      copy_from(other);
    }

    void operator= (const arr &&other)
    {
      copy_from(std::move(other));
    }

    void operator& () = delete;

    auto get() & // lvalue
    {
      std::array<atype,N> b;
      for (int i=0; i<N; i++) b[i] = elem[i].get();
      return b;
    }

    auto get() && // rvalue
    {
      std::array<atype,N> b;
      for (u64 i=0; i<N; i++) b[i] = std::move(elem[i]).get();
      return b;
    }

    void set_time(u64 t)
    {
      for (u64 i=0; i<N; i++) elem[i].set_time(t);
    }

  public:

    u64 time() const
    {
      u64 t = 0;
      for (u64 i=0; i<N; i++) {
	t = std::max(t,elem[i].time());
      }
      return t;
    }

    arr() {}

    arr(arr &other)
    {
      copy_from(other);
    }

    arr(const arr &&other)
    {
      copy_from(std::move(other));
    }

    template<std::convertible_to<T> ...U>
    arr(U&&... args) : elem{std::forward<U>(args)...} {}

    arr(unaryfunc<u64,T> auto f)
    {
      for (u64 i=0; i<N; i++) {
	elem[i] = f(i);
      }
    }

    template<std::convertible_to<T> U>
    arr(std::array<U,N> &b)
    {
      for (u64 i=0; i<N; i++) {
	elem[i] = b[i];
      }
    }

    template<std::convertible_to<T> U>
    arr(U (&b)[N])
    {
      for (u64 i=0; i<N; i++) {
	elem[i] = b[i];
      }
    }

    template<arrtype U>
    arr(U &&other)
    {
      copy_from(std::forward<U>(other));
    }

    void operator= (arr &other) requires (regtype<T>)
    {
      copy_from(other);
    }

    void operator= (const arr &&other) requires (regtype<T>)
    {
      copy_from(std::move(other));
    }

    template<arrtype U> requires (regtype<T>)
    void operator= (U &&other)
    {
      copy_from(std::forward<U>(other));
    }

    T& operator[] (u64 i)
    {
      assert(i<N);
      return elem[i];
    }

    template<valtype U>
    [[nodiscard]] valt<T> select(U && index)
    {
      // only for reading
      // we do not bother providing an rvalue version
      static_assert(ival<U>);
      static_assert(N>=1);
      if constexpr (N>=2) {
	constexpr auto c = MUX<N,T::size>;
	panel.update_logic(c[0]); // MUX select
	panel.update_logic(c[1]); // MUX data
	auto [i,ti] = std::forward<U>(index).get_vt();
	assert(i<N);
	auto [d,td] = elem[i].get_vt(); // lvalue
	auto t = std::max(ti+c[0].delay(),td) + c[1].delay();
	return {d,t};
      } else {
	static_assert(N==1);
	return elem[0];
      }
    }

    template<arith auto FO>
    void fanout(hard<FO>) & // lvalue
    {
      static_assert(FO>=2);
      for (auto &e : elem) e.fanout(hard<FO>{});
    }

    [[nodiscard]] arr&& fo1() & // lvalue
    {
      return std::move(*this);
    }

    void fanout(auto) && = delete; // rvalue, not needed

    void print(std::string before="", std::string after="\n", bool t=true, std::ostream & os=std::cout) const
    {
      for (u64 i=0; i<N; i++) elem[i].print(before+std::to_string(i)+": ",after,t,os);
    }

    void printb(std::string before="", std::string after="\n", bool t=true, std::ostream & os=std::cout) const
    {
      for (u64 i=0; i<N; i++) elem[i].printb(before+std::to_string(i)+": ",after,t,os);
    }

    auto concat() & requires std::unsigned_integral<atype> // lvalue
    {
      // element 0 is at rightmost position
      static_assert(N!=0);
      u64 y = 0;
      for (i64 i=N-1; i>=0; i--) {
	y = (y<<T::size) | elem[i].get();
      }
      return val<N*T::size> {y,time()};
    }

    auto concat() && requires std::unsigned_integral<atype> // rvalue
    {
      // element 0 is at rightmost position
      static_assert(N!=0);
      u64 y = 0;
      for (i64 i=N-1; i>=0; i--) {
	y = (y<<T::size) | std::move(elem[i]).get();
      }
      return val<N*T::size> {y,time()};
    }

    template<std::convertible_to<valt<T>> U>
    [[nodiscard]] auto append(U &&x) & // lvalue
    {
      arr<valt<T>,N+1> b;
      for (u64 i=0; i<N; i++) {
	b[i] = elem[i];
      }
      b[N] = std::forward<U>(x);
      return b;
    }

    template<std::convertible_to<valt<T>> U>
    [[nodiscard]] auto append(U &&x) && // rvalue
    {
      arr<valt<T>,N+1> b;
      for (u64 i=0; i<N; i++) {
	b[i] = std::move(elem[i]);
      }
      b[N] = std::forward<U>(x);
      return b;
    }

    template<u64 W>
    [[nodiscard]] auto make_array(val<W>&&) & // lvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(W!=0 && W<=64);
      constexpr u64 NBITS = T::size*N;
      constexpr u64 M = (NBITS+W-1)/W;
      auto data = get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      auto aa = unpack_bits<W>(a);
      static_assert(aa.size()>=M);
      arr<val<W>,M> out {aa};
      out.set_time(t);
      return out;
    }

    template<u64 W>
    [[nodiscard]] auto make_array(val<W>&&) && // rvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(W!=0 && W<=64);
      constexpr u64 NBITS = T::size*N;
      constexpr u64 M = (NBITS+W-1)/W;
      auto data = std::move(*this).get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      auto aa = unpack_bits<W>(a);
      static_assert(aa.size()>=M);
      arr<val<W>,M> out {aa};
      out.set_time(t);
      return out;
    }

    template<valtype U>
    [[nodiscard]] auto shift_left(U && x) & // lvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(std::unsigned_integral<base<U>>);
      static_assert(N!=0);
      static_assert(U::size!=0);
      auto data = get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      if constexpr (U::size == 64) {
	for (u64 i=a.size(); i!=0; i--) a[i] = a[i]-1;
	a[0] = std::forward<U>(x).get();
      } else {
	static_assert(U::size<64);
	for (u64 i=a.size()-1; i!=0; i--) {
	  a[i] = (a[i] << U::size) | (a[i-1] >> (64-U::size));
	}
	a[0] = (a[0] << U::size) | (std::forward<U>(x).get() & ((u64(1)<<U::size)-1));
      }
      auto aa = unpack_bits<T::size>(a);
      static_assert(aa.size()>=N);
      arr<valt<T>,N> out {aa};
      out.set_time(t);
      return out;
    }

    template<valtype U>
    [[nodiscard]] auto shift_left(U && x) && // rvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(std::unsigned_integral<base<U>>);
      static_assert(N!=0);
      static_assert(U::size!=0);
      auto data = std::move(*this).get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      if constexpr (U::size == 64) {
	for (u64 i=a.size(); i!=0; i--) a[i] = a[i]-1;
	a[0] = std::forward<U>(x).get();
      } else {
	static_assert(U::size<64);
	for (u64 i=a.size()-1; i!=0; i--) {
	  a[i] = (a[i] << U::size) | (a[i-1] >> (64-U::size));
	}
	a[0] = (a[0] << U::size) | (std::forward<U>(x).get() & ((u64(1)<<U::size)-1));
      }
      auto aa = unpack_bits<T::size>(a);
      static_assert(aa.size()>=N);
      arr<valt<T>,N> out {aa};
      out.set_time(t);
      return out;
    }

    template<valtype U>
    [[nodiscard]] auto shift_right(U && x) & // lvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(std::unsigned_integral<base<U>>);
      static_assert(N!=0);
      static_assert(U::size!=0);
      auto data = get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      if constexpr (U::size == 64) {
	for (u64 i=0; i<a.size()-1; i++) a[i] = a[i]+1;
	a[a.size()-1] = 0;
      } else {
	static_assert(U::size<64);
	for (u64 i=0; i<a.size()-1; i++) {
	  a[i] = (a[i] >> U::size) | (a[i+1] << (64-U::size));
	}
	a[a.size()-1] >>= U::size;
      }
      auto xdata = std::forward<U>(x).get();
      u64 shiftin = (U::size==64)? xdata : xdata & ((u64(1)<<U::size)-1);
      u64 k = T::size*N - U::size;
      u64 j = k/64;
      u64 pos = k%64;
      a[j] |= shiftin << pos;
      if (pos+U::size > 64) {
	a[j+1] = shiftin >> (64-pos);
      }
      auto aa = unpack_bits<T::size>(a);
      static_assert(aa.size()>=N);
      arr<valt<T>,N> out {aa};
      out.set_time(t);
      return out;
    }

    template<valtype U>
    [[nodiscard]] auto shift_right(U && x) && // rvalue
    {
      static_assert(std::unsigned_integral<atype>);
      static_assert(std::unsigned_integral<base<U>>);
      static_assert(N!=0);
      static_assert(U::size!=0);
      auto data = std::move(*this).get();
      auto t = time(); // FIXME?
      auto a = pack_bits<T::size>(data);
      if constexpr (U::size == 64) {
	for (u64 i=0; i<a.size()-1; i++) a[i] = a[i]+1;
	a[a.size()-1] = 0;
      } else {
	static_assert(U::size<64);
	for (u64 i=0; i<a.size()-1; i++) {
	  a[i] = (a[i] >> U::size) | (a[i+1] << (64-U::size));
	}
	a[a.size()-1] >>= U::size;
      }
      auto xdata = std::forward<U>(x).get();
      u64 shiftin = (U::size==64)? xdata : xdata & ((u64(1)<<U::size)-1);
      u64 k = T::size*N - U::size;
      u64 j = k/64;
      u64 pos = k%64;
      a[j] |= shiftin << pos;
      if (pos+U::size > 64) {
	a[j+1] = shiftin >> (64-pos);
      }
      auto aa = unpack_bits<T::size>(a);
      static_assert(aa.size()>=N);
      arr<valt<T>,N> out {aa};
      out.set_time(t);
      return out;
    }

    valt<T> fold_xor() & // lvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = XOR<N> * T::size;
	panel.update_logic(c);
	auto data = get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x ^= e;
	return {x,t};
      } else if constexpr (N==1) {
	return elem[0];
      } else {
	return 0;
      }
    }

    valt<T> fold_xor() && // rvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = XOR<N> * T::size;
	panel.update_logic(c);
	auto data = std::move(*this).get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x ^= e;
	return {x,t};
      } else if constexpr (N==1) {
	return std::move(elem[0]);
      } else {
	return 0;
      }
    }

    valt<T> fold_or() & // lvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = OR<N> * T::size;
	panel.update_logic(c);
	auto data = get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x |= e;
	return {x,t};
      } else if constexpr (N==1) {
	return elem[0];
      } else {
	return 0;
      }
    }

    valt<T> fold_or() && // rvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = OR<N> * T::size;
	panel.update_logic(c);
	auto data = std::move(*this).get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x |= e;
	return {x,t};
      } else if constexpr (N==1) {
	return std::move(elem[0]);
      } else {
	return 0;
      }
    }    

    valt<T> fold_and() & // lvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = AND<N> * T::size;
	panel.update_logic(c);
	auto data = get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x &= e;
	return {x,t};
      } else if constexpr (N==1) {
	return elem[0];
      } else {
	return 0;
      }
    }

    valt<T> fold_and() && // rvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr circuit c = AND<N> * T::size;
	panel.update_logic(c);
	auto data = std::move(*this).get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x &= e;
	return {x,t};
      } else if constexpr (N==1) {
	return std::move(elem[0]);
      } else {
	return 0;
      }
    }

    auto fold_add() & // lvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr u64 RBITS = T::size + std::bit_width(N-1); // output bits
	constexpr circuit c = ADDN<T::size,N>;
	panel.update_logic(c);
	auto data = get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x += e;
	return val<RBITS,atype>{x,t};
      } else if constexpr (N==1) {
	return valt<T>{elem[0]};
      } else {
	return valt<T>{0};
      }
    }

    auto fold_add() && // rvalue
    {
      static_assert(ival<T>);
      if constexpr (N>=2) {
	constexpr u64 RBITS = T::size + std::bit_width(N-1); // output bits
	constexpr circuit c = ADDN<T::size,N>;
	panel.update_logic(c);
	auto data = std::move(*this).get();
	auto t = time() + c.delay();
	atype x = 0;
	for (auto e : data) x += e;
	return val<RBITS,atype>{x,t};
      } else if constexpr (N==1) {
	return valt<T>{elem[0]};
      } else {
	return valt<T>{0};
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
  public:
    using type = T;
    using valuetype = rawdata<T>::type;
    using static_ram = sram<N,rawdata<T>::width>;

  private:
    valuetype data[N];
    u64 last_read_cycle = 0;
    u64 last_write_cycle = 0;

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

    ram(ram &) = delete;
    ram& operator= (ram&) = delete;
    void operator& () = delete;

  public:

    ram()
    {
      panel.update_storage(static_ram::NBITS,static_ram::XTORS);
    }

    ~ram()
    {
      panel.storage_destroyed = true;
      while (! writes.empty()) {
	writes[0].commit(*this);
	std::pop_heap(writes.begin(),writes.end());
	writes.pop_back();
      }
    }

    template<ival TYPEA, std::convertible_to<T> TYPED>
    void write(TYPEA && address, TYPED && dataval)
    {
      assert(panel.cycle>last_write_cycle && "single RAM write per cycle");
      last_write_cycle = panel.cycle;
      if (! exec.active) return;
      panel.update_energy(static_ram::EWRITE);
      auto [av,at] = proxy::get_vt(std::forward<TYPEA>(address));
      auto dv = proxy::get(std::forward<TYPED>(dataval));
      auto dt = proxy::time(dataval);
      writes.push_back({u64(av),valuetype(dv),std::max(at,dt)});
      std::push_heap(writes.begin(),writes.end());
    }

    template<ival U>
    T read(U && address)
    {
      assert(panel.cycle>last_read_cycle && "single RAM read per cycle");
      last_read_cycle = panel.cycle;
      if (! exec.active) return {};
      panel.update_energy(static_ram::EREAD);
      auto [addr,t] = proxy::get_vt(std::forward<U>(address));
      while (! writes.empty() && writes[0].t <= t) {
	writes[0].commit(*this);
	std::pop_heap(writes.begin(),writes.end());
	writes.pop_back();
      }
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
    val<M,T> a[N]; // TODO: make it raw data instead
    
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

    template<ival U>
    val<M,T> read(U && address)
    {
      auto [i,t] = proxy::get_vt(std::forward<U>(address));
      assert(i<N);
      //TODO: energy and delay
      return {proxy::get(a[i]),t};
    }    
  };


  // ###########################

  // EQUALITY
  template<valtype T1, valtype T2>
  val<1> operator== (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = EQUAL<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_equal(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator== (T1&& x1, T2 x2)
  {
    constexpr circuit reduc = NOR<valt<T1>::size>;
    const circuit c = INV * ones<valt<T1>::size>(x2) + reduc; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_equal(v1,x2), t1+c.delay()};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator== (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) == x1;
  }

  // INEQUALITY
  template<valtype T1, valtype T2>
  val<1> operator!= (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = NEQ<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_different(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, arith T2> // second argument is a constant
  val<1> operator!= (T1&& x1, T2 x2)
  {
    constexpr circuit reduc = OR<valt<T1>::size>;
    const circuit c = INV * ones<valt<T1>::size>(x2) + reduc; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_different(v1,x2), t1+c.delay()};
  }

  template<arith T1, valtype T2> // first argument is a constant
  val<1> operator!= (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) != x1;
  }

  // GREATER THAN
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  val<1> operator> (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = GT<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_greater(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // second argument is a constant
  val<1> operator> (T1&& x1, T2 x2)
  {
    constexpr u64 N = valt<T1>::size;
    static_assert(N!=0);
    constexpr circuit comp = GT<N>;
    constexpr circuit comp0 = (N==1)? circuit{} : (std::signed_integral<base<T1>>)? NOR<N-1> + NOR<2> : OR<N>;
    const circuit &c = (x2==0)? comp0 : comp; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_greater(v1,x2), t1+c.delay()};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // first argument is a constant
  val<1> operator> (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) < x1;
  }

  // LESS THAN
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  val<1> operator< (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = GT<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_less(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // second argument is a constant
  val<1> operator< (T1&& x1, T2 x2)
  {
    constexpr circuit comp = GT<valt<T1>::size>;
    constexpr circuit comp0;
    const circuit &c = (x2==0)? comp0 : comp; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_less(v1,x2), t1+c.delay()};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // first argument is a constant
  val<1> operator< (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) > x1;
  }

  // GREATER THAN OR EQUAL
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  val<1> operator>= (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = GTE<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_greater_equal(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // second argument is a constant
  val<1> operator>= (T1&& x1, T2 x2)
  {
    constexpr circuit comp = GTE<valt<T1>::size>;
    constexpr circuit comp0 = (std::signed_integral<base<T1>>)? INV : circuit{};
    const circuit &c = (x2==0)? comp0 : comp; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_greater_equal(v1,x2), t1+c.delay()};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // first argument is a constant
  val<1> operator>= (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) <= x1;
  }
  
  // LESS THAN OR EQUAL
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  val<1> operator<= (T1&& x1, T2&& x2)
  {
    static_assert(valt<T1>::size == valt<T2>::size);
    constexpr circuit c = GTE<valt<T1>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {is_less_equal(v1,v2), std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // second argument is a constant
  val<1> operator<= (T1&& x1, T2 x2)
  {
    constexpr u64 N = valt<T1>::size;
    static_assert(N!=0);
    constexpr circuit comp = GTE<N>;
    constexpr circuit comp0 = (N==1)? circuit{} : (std::signed_integral<base<T1>>)? NOR<N-1> + OR<2> : NOR<N>;
    const circuit &c = (x2==0)? comp0 : comp; // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {is_less_equal(v1,x2), t1+c.delay()};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // first argument is a constant
  val<1> operator<= (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) >= x1;
  }

  // ADDITION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  auto operator+ (T1&& x1, T2&& x2)
  {
    constexpr circuit c = ADD<valt<T1,T2>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    using rtype = val<std::min(val<64>::size,valt<T1,T2>::size+1),decltype(v1+v2)>;
    return rtype{v1+v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  auto operator+ (T1&& x1, T2 x2)
  {
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    using rtype = val<std::min(val<64>::size,valt<T1,T2>::size+1),decltype(v1+x2)>;
    if (x2==0) return rtype{v1,t1};
    constexpr circuit c = INC<valt<T1>::size>;
    proxy::update_logic(c);
    return rtype{v1+x2, t1+c.delay()};
  }  

  template<std::integral T1, valtype T2> requires (ival<T2>) // 1st arg constant 
  auto operator+ (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) + x1;
  }

  // CHANGE SIGN
  template<valtype T>
  auto operator- (T&& x)
  {
    constexpr circuit c = INV * valt<T>::size + INC<valt<T>::size>;
    proxy::update_logic(c);
    auto [v,t] = proxy::get_vt(std::forward<T>(x));
    using rtype = val<valt<T>::size,decltype(-v)>;
    return rtype{-v, t+c.delay()};
  }

  // SUBTRACTION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  auto operator- (T1&& x1, T2&& x2)
  {
    constexpr circuit c = SUB<valt<T1,T2>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    using rtype = val<std::min(val<64>::size,valt<T1,T2>::size+1),decltype(v1-v2)>;
    return rtype{v1-v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  auto operator- (T1&& x1, T2 x2)
  {
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    using rtype = val<std::min(val<64>::size,valt<T1,T2>::size+1),decltype(v1-x2)>;
    if (x2==0) return rtype{v1,t1};
    constexpr circuit c = INC<valt<T1>::size>;
    proxy::update_logic(c);
    return rtype{v1-x2, t1+c.delay()};
  }

  template<std::integral T1, valtype T2> requires (ival<T2>) // 1st arg constant 
  auto operator- (T1 x1, T2&& x2)
  {
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    using rtype = val<std::min(val<64>::size,valt<T1,T2>::size+1),decltype(x1-v2)>;
    if (x1==0) return rtype{-v2,t2};
    constexpr circuit c = INC<valt<T2>::size>;
    proxy::update_logic(c);
    return rtype{x1-v2, t2+c.delay()};
  }

  // LEFT SHIFT
  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  valt<T1> operator<< (T1&& x1, T2 x2)
  {
    // no transistors
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1<<x2, t1};
  }

  // RIGHT SHIFT
  template<valtype T1, std::integral T2> requires (ival<T1>) // 2nd arg constant
  valt<T1> operator>> (T1&& x1, T2 x2)
  {
    // no transistors
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1>>x2, t1};
  }

  // MULTIPLICATION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  auto operator* (T1&& x1, T2&& x2)
  {
    constexpr circuit c = IMUL<valt<T1>::size,valt<T2>::size>;
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    using rtype = val<valt<T1>::size+valt<T2>::size,decltype(v1*v2)>;
    return rtype{v1*v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, arith auto N2> requires (ival<T1>) // 2nd argument is hard value
  auto operator* (T1&& x1, hard<N2>)
  {
    // multiply first argument by fixed, known multiplier N
    constexpr u64 N = (N2>=0)? N2 : truncate<std::bit_width(u64(-N2-1))+1>(N2);
    constexpr circuit c = HIMUL<N,valt<T1>::size>; // FIXME: signed multiplication
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    using rtype = val<valt<T1>::size+std::bit_width(N),decltype(v1*N2)>;
    return rtype{v1*N2, t1+c.delay()};
  }

  template<arith auto N1, valtype T2> requires (ival<T2>) // 1st argument is hard value
  auto operator* (hard<N1> x1, T2&& x2)
  {
    return x2 * x1;
  }
  
  // DIVISION
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T1> operator/ (T1&& x1, T2&& x2)
  {
    // TODO
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {v1/v2, std::max(t1,t2)};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // divisor is a constant
  valt<T1> operator/ (T1&& x1, T2 x2)
  {
    // TODO
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1/x2, t1};
  }
  
  // MODULO (REMAINDER)
  template<valtype T1, valtype T2> requires (ival<T1> && ival<T2>)
  valt<T2> operator% (T1&& x1, T2&& x2)
  {
    // TODO
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {v1%v2, std::max(t1,t2)};
  }

  template<valtype T1, std::integral T2> requires (ival<T1>) // modulus is a constant
  valt<T2> operator% (T1&& x1, T2 x2)
  {
    // TODO
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1%x2, t1};
  }

  // BITWISE AND
  template<valtype T1, valtype T2>
  valt<T1,T2> operator& (T1&& x1, T2&& x2)
  {
    constexpr circuit c = AND<2> * std::min(valt<T1>::size,valt<T2>::size);
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {v1&v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> // second argument is constant
  valt<T1> operator& (T1&& x1, T2 x2)
  {
    // no transistors
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1&x2, t1};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator& (T1 x1, T2&& x2)
  {
    // no transistors
    return std::forward<T2>(x2) & x1;
  }
  
  // BITWISE OR
  template<valtype T1, valtype T2>
  valt<T1,T2> operator| (T1&& x1, T2&& x2)
  {
    constexpr circuit c = OR<2> * std::min(valt<T1>::size,valt<T2>::size);
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {v1|v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> // second argument is constant
  valt<T1> operator| (T1&& x1, T2 x2)
  {
    // no transistors
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1|x2, t1};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator| (T1 x1, T2&& x2)
  {
    // no transistors
    return std::forward<T2>(x2) | x1;
  }

  // BITWISE EXCLUSIVE OR
  template<valtype T1, valtype T2>
  valt<T1,T2> operator^ (T1&& x1, T2&& x2)
  {
    constexpr circuit c = XOR<2> * std::min(valt<T1>::size,valt<T2>::size);
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    return {v1^v2, std::max(t1,t2)+c.delay()};
  }

  template<valtype T1, std::integral T2> // second argument is constant
  valt<T1> operator^ (T1&& x1, T2 x2)
  {
    const circuit c = INV * ones<valt<T1>::size>(x2); // not constexpr
    proxy::update_logic(c);
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    return {v1^x2, t1+c.delay()};
  }

  template<std::integral T1, valtype T2> // first argument is constant
  valt<T2> operator^ (T1 x1, T2&& x2)
  {
    return std::forward<T2>(x2) ^ x1;
  }

  // BITWISE COMPLEMENT
  template<valtype T>
  valt<T> operator~ (T&& x)
  {
    constexpr circuit c = INV * valt<T>::size;
    proxy::update_logic(c);
    auto [v,t] = proxy::get_vt(std::forward<T>(x));
    return {~v, t+c.delay()};
  }

  // CONCATENATE BITS
  template<valtype T1, valtype T2, valtype... T>
  auto concat(T1&& x1, T2&& x2, T&&... x)
  {
    // no transistors (wires not modeled, TODO?)
    if constexpr (sizeof...(x)==0) {
      return proxy::concat2(std::forward<T1>(x1),std::forward<T2>(x2));
    } else {
      return proxy::concat2(proxy::concat2(std::forward<T1>(x1),std::forward<T2>(x2)),std::forward<T>(x)...);
    }
  }

  // SELECT BETWEEN TWO VALUES
  template<valtype T, valtype T1, valtype T2>
  valt<T1,T2> select(T &&cond, T1 &&x1, T2 &&x2)
  {
    // this is NOT conditional execution: both sides are evaluated
    static_assert(valt<T>::size == 1,"the condition of a select is a single bit");
    static_assert(valt<T1>::size == valt<T2>::size,"both sides of a select must have the same size");
    constexpr auto c = MUX<2,valt<T1>::size>;
    proxy::update_logic(c[0]); // MUX select
    proxy::update_logic(c[1]); // MUX data
    auto [vc,tc] = proxy::get_vt(std::forward<T>(cond));
    auto [v1,t1] = proxy::get_vt(std::forward<T1>(x1));
    auto [v2,t2] = proxy::get_vt(std::forward<T2>(x2));
    auto t = std::max({tc+c[0].delay(),t1,t2}) + c[1].delay();
    if (bool(vc)) {
      return {v1,t};
    } else {
      return {v2,t};
    }
  }

  // CONDITIONAL EXECUTION
  template<valtype T, action A> requires (std::same_as<return_type<A>,void>)
  void execute(T &&mask, const A &f)
  {
    static_assert(std::unsigned_integral<base<T>>);
    auto prev_exec = exec;
    auto [m,t] = proxy::get_vt(std::forward<T>(mask));
    for (u64 i=0; i<valt<T>::size; i++) {
      bool cond = (m>>i) & 1;
      exec.set_state(cond,t);
      // to prevent cheating, we execute the action even when the condition is false
      // (otherwise, this primitive could be used to leak any bit)
      if constexpr (std::invocable<A>) {
	f();
      } else {
	static_assert(std::invocable<A,u64>);
	f(i);
      }
    }
    exec = prev_exec;
  }


  template<valtype T, action A>
  auto execute(T &&mask, const A &f)
  {
    static_assert(valtype<return_type<A>>);
    static_assert(std::unsigned_integral<base<T>>);
    constexpr u64 N = valt<T>::size;
    using rtype = valt<return_type<A>>;
    auto prev_exec = exec;
    auto [m,t] = proxy::get_vt(std::forward<T>(mask));
    arr<rtype,N> result;
    for (u64 i=0; i<N; i++) {
      bool cond = (m>>i) & 1;
      exec.set_state(cond,t);
      // to prevent cheating, we execute the action even when the condition is false
      // (otherwise, this primitive could be used to leak any bit)
      if constexpr (std::invocable<A>) {
	result[i] = f();
      } else {
	static_assert(std::invocable<A,u64>);
	result[i] = f(i);
      }
    };
    exec = prev_exec;
    return result;
  }

  // ###########################
}

#endif // HARCOM_H

