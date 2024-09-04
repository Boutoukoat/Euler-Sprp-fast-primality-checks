// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// fastest code for 65-bits modulus
//
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "m64_utils.h"
#include "optimized.h"
#include "m_reg.h"

typedef unsigned __int128 uint128_t;

// subtract the modulus 'mod' multiple times from the input number 'res', if needed
static inline void ciosSubtract(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo)
{
	uint64_t n_lo, n_hi;
	uint64_t t_lo, t_hi;
	uint8_t b;
	n_lo = *res_lo;
	n_hi = *res_hi;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		t_hi = n_hi;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
		b = my_sbb64(b, n_hi, 1, &n_hi);
	}
	while (b == 0);
	// get the saved values after a borrow occured
	*res_lo = t_lo;
	*res_hi = t_hi;
}

// subtract the modulus 'mod' multiple times from the input number 'res', if needed
static inline void ciosSubtract(uint64_t * res_lo, uint64_t * res_hi, uint64_t carries, uint64_t mod_lo)
{
	uint64_t n_lo, n_hi;
	uint64_t t_lo, t_hi;
	uint8_t b;
	n_lo = *res_lo;
	n_hi = *res_hi;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		t_hi = n_hi;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
		b = my_sbb64(b, n_hi, 1, &n_hi);
		if (__builtin_constant_p(carries) && carries == 0) {
		} else {
			b = my_sbb64(b, carries, 1, &carries);
		}
	}
	while (b == 0);
	// get the saved values after a borrow occured
	*res_lo = t_lo;
	*res_hi = t_hi;
}

static void ciosConstants(uint64_t x, uint64_t mod_lo, uint64_t * one_lo, uint64_t * one_hi, uint64_t * res_lo, uint64_t * res_hi)
{
	uint128_t mod = ((uint128_t) 1 << 64) + mod_lo;
#if PARANOID
	assert(x < 64);
#endif

#if 1
	uint64_t t0, t1, t2, s0, s1, s2, s;
	uint8_t c;

	uint64_t barrett = my_div64((uint64_t) - 1 >> 1, (uint64_t) - 1, (uint64_t) (mod >> 1), &s);

	// t = barrett * mod
	t0 = my_mul64(barrett, mod_lo, &t1);
	t2 = my_adc64(0, t1, barrett, &t1);
	if (t2) {
		// t = barrett * mod
		barrett -= 1;
		c = my_sbb64(0, t0, mod_lo, &t0);
		c = my_sbb64(c, t1, 1, &t1);
		my_sbb64(c, t2, 0, &t2);
	}
#if PARANOID
	assert(t2 == 0);
#endif

	// computes (1 << 128) %  mod
	// 2^128 - t
	c = my_sbb64(0, 0, t0, &t0);
	my_sbb64(c, 0, t1, &t1);
	// one = 2^128 % m
	*one_lo = t0;
	*one_hi = t1;

#if PARANOID
	// check the number is less than the modulus
	if (t1 == 1)
		assert(t0 < mod_lo);
#endif

	// computes ((1<<x) << 128) %  mod

	// temp = one << x
	if (!__builtin_constant_p(x) || (__builtin_constant_p(x) && x != 0)) {
		my_shld64(&t2, &t1, &t0, x);
	}

	if (__builtin_constant_p(x) && x < 7) {
		// no reduction, number is small enough
	} else {
		// s = barrett * temp_hi 
		s0 = my_mul64(barrett, t1, &s1);
		s = my_mul64(barrett, t2, &s2);
		c = my_adc64(0, s1, s, &s1);
		// s = s * modulus
		s = s1;
		s0 = my_mul64(mod_lo, s, &s1);
		s2 = my_adc64(0, s1, s, &s1);
		// temp = temp - s      
		c = my_sbb64(0, t0, s0, &t0);
		c = my_sbb64(c, t1, s1, &t1);
	}

	// res = (1 << x)*(2^128) % m
	*res_lo = t0;
	*res_hi = t1;
#else
	uint128_t t;

	// computes (1 << 128) %  mod
	t = (uint128_t) 0;	// 2^128       (ignore the most significant bit = 1)
	t -= mod;		// 2^128-m     (ignore the borrow = 1) 
	t %= mod;		// (2^128-m) % m     is 2^128 % m
	*one_lo = (uint64_t) t;
	*one_hi = (uint64_t) (t >> 64);

	// computes ((1<<x) << 128) %  mod
	t <<= x;
	if ((__builtin_constant_p(x) && x > 7) || (x > 7)) {
		t %= mod;
	}
	*res_lo = (uint64_t) t;
	*res_hi = (uint64_t) (t >> 64);

#endif
}

static inline __attribute__((always_inline))
void ciosModSquare(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mmagic)
{

#if INLINE_ASM && defined(__x86_64__) && defined(__BMI2__)
	uint64_t n_lo = *res_lo, n_hi = *res_hi;
	uint64_t t0, t1, t2, cc, cs, ct;
 asm("\
    movq    %[n_lo], %%rdx              \n\
    mulxq   %[n_lo], %[t0], %[cc]       \n\
    leaq    (%[n_hi],%[n_hi]), %%rdx    \n\
    mulxq   %[n_lo], %[t1], %[t2]       \n\
    movq    %[mmagic], %%rdx            \n\
    imulq   %[t0], %%rdx                \n\
    mulxq   %[mod_lo], %[cs], %[ct]     \n\
    addq    %[cc], %[t1]                \n\
    adcq    $0, %[t2]                   \n\
    addq    %[t0], %[cs]                \n\
    adcq    %%rdx, %[ct]                \n\
    setb    %b[t0]                      \n\
    movzb   %b[t0], %k[t0]              \n\
    movq    %[n_hi], %%rdx              \n\
    mulxq   %[n_hi], %[cs], %[cc]       \n\
    addq    %[t1], %[ct]                \n\
    adcq    %[t2], %[t0]                \n\
    movq    %[ct], %%rdx                \n\
    imulq   %[mmagic], %%rdx            \n\
    mulxq   %[mod_lo], %[n_hi], %[n_lo] \n\
    addq    %[t0], %[cs]                \n\
    adcq    $0, %[cc]                   \n\
    addq    %[ct], %[n_hi]              \n\
    adcq    %%rdx, %[n_lo]              \n\
    setb    %b[n_hi]                    \n\
    movzb   %b[n_hi], %k[n_hi]          \n\
    addq    %[cs], %[n_lo]              \n\
    adcq    %[cc], %[n_hi]              \n":"=&r"(n_lo), "=&r"(n_hi),[t0] "=&r"(t0),[t1] "=&r"(t1),[t2] "=&r"(t2),[cc] "=&r"(cc),[cs] "=&r"(cs),[ct] "=&r"(ct)
 :	    [n_lo] "0"(n_lo),[n_hi] "1"(n_hi),[mmagic] "r"(mmagic),[mod_lo] "r"(mod_lo)
 :	    "flags", "%rdx");

	*res_lo = n_lo;
	*res_hi = n_hi;

#else

	uint64_t n_lo = *res_lo, n_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, m;

	cc = (uint128_t) n_lo *n_lo;	// #1
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	cc += (uint128_t) n_lo *(n_hi + n_hi);	// #2
	t1 = (uint64_t) cc;
	cc = cc >> 64;
	t2 = (uint64_t) cc;

	m = t0 * mmagic;	// #3
	cs = (uint128_t) m *mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += m;
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t) cs;

	cc = (uint128_t) n_hi *n_hi;	// #5
	cc += t1;
	t1 = (uint64_t) cc;
	cc = cc >> 64;
	t2 = (uint64_t) cc;

	m = t0 * mmagic;	// #6
	cs = (uint128_t) m *mod_lo;	// #7
	cs += t0;
	cs = cs >> 64;
	cs += m;
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t) cs;

	*res_lo = t0;
	*res_hi = t1;
#endif
}

static inline __attribute__((always_inline))
void ciosModSquare3(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mmagic)
{
#if INLINE_ASM && defined(__x86_64__) && defined(__BMI2__)
	uint64_t n_lo = *res_lo, n_hi = *res_hi;
	uint64_t t0, t1, t2, t3, t4, t5, t6;
 asm("\
    movq    %[n_lo], %%rdx                \n\
    mulxq   %[n_lo], %[t0], %[t1]         \n\
    leaq    (%[n_hi],%[n_hi]), %%rdx      \n\
    mulxq   %[n_lo], %[n_lo], %[t3]       \n\
    leaq    (%[n_lo],%[t1]), %[t4]        \n\
    movq    %[mmagic], %[t6]              \n\
    imulq   %[t0], %[t6]                  \n\
    movq    %[t6], %%rdx                  \n\
    mulxq   %[mod_lo], %%rdx, %[t5]       \n\
    xorl    %k[t2], %k[t2]                \n\
    addq    %[t4], %[t6]                  \n\
    setb    %b[t2]                        \n\
    addq    %[t0], %%rdx                  \n\
    adcq    %[t5], %[t6]                  \n\
    adcq    $0, %[t2]                     \n\
    addq    %[t1], %[n_lo]                \n\
    adcq    %[t3], %[t2]                  \n\
    movq    %[n_hi], %%rdx                \n\
    mulxq   %[n_hi], %[n_lo], %[t4]       \n\
    addq    %[t2], %[n_lo]                \n\
    adcq    $0, %[t4]                     \n\
    movq    %[mmagic], %%rdx              \n\
    imulq   %[t6], %%rdx                  \n\
    mulxq   %[mod_lo], %[t0], %[t1]       \n\
    xorl    %k[n_hi], %k[n_hi]            \n\
    addq    %[n_lo], %%rdx                \n\
    setb    %b[n_hi]                      \n\
    addq    %[t6], %[t0]                  \n\
    adcq    %[t1], %%rdx                  \n\
    mulxq   %%rdx, %[t1], %[t5]           \n\
    adcq    %[t4], %[n_hi]                \n\
    leaq    (%[n_hi],%[n_hi]), %[t6]      \n\
    mulxq   %[t6], %[t0], %[t3]           \n\
    leaq    (%[t0],%[t5]), %[t4]          \n\
    movq    %[mmagic], %[t6]              \n\
    imulq   %[t1], %[t6]                  \n\
    movq    %[t6], %%rdx                  \n\
    mulxq   %[mod_lo], %%rdx, %[t2]       \n\
    xorl    %k[n_lo], %k[n_lo]            \n\
    addq    %[t4], %[t6]                  \n\
    setb    %b[n_lo]                      \n\
    addq    %[t1], %%rdx                  \n\
    adcq    %[t2], %[t6]                  \n\
    adcq    $0, %[n_lo]                   \n\
    addq    %[t5], %[t0]                  \n\
    adcq    %[t3], %[n_lo]                \n\
    movq    %[n_hi], %%rdx                \n\
    mulxq   %[n_hi], %[t4], %[t1]         \n\
    addq    %[n_lo], %[t4]                \n\
    adcq    $0, %[t1]                     \n\
    movq    %[mmagic], %%rdx              \n\
    imulq   %[t6], %%rdx                  \n\
    mulxq   %[mod_lo], %[t0], %[n_lo]     \n\
    xorl    %k[n_hi], %k[n_hi]            \n\
    addq    %[t4], %%rdx                  \n\
    setb    %b[n_hi]                      \n\
    addq    %[t6], %[t0]                  \n\
    adcq    %[n_lo], %%rdx                \n\
    mulxq   %%rdx, %[t4], %[t5]           \n\
    adcq    %[t1], %[n_hi]                \n\
    leaq    (%[n_hi],%[n_hi]), %[t6]      \n\
    mulxq   %[t6], %[t0], %[t3]           \n\
    leaq    (%[t0],%[t5]), %[t1]          \n\
    movq    %[mmagic], %[t6]              \n\
    imulq   %[t4], %[t6]                  \n\
    movq    %[t6], %%rdx                  \n\
    mulxq   %[mod_lo], %%rdx, %[t2]       \n\
    xorl    %k[n_lo], %k[n_lo]            \n\
    addq    %[t1], %[t6]                  \n\
    setb    %b[n_lo]                      \n\
    addq    %[t4], %%rdx                  \n\
    adcq    %[t2], %[t6]                  \n\
    adcq    $0, %[n_lo]                   \n\
    addq    %[t5], %[t0]                  \n\
    adcq    %[t3], %[n_lo]                \n\
    movq    %[n_hi], %%rdx                \n\
    mulxq   %[n_hi], %[t4], %[t1]         \n\
    addq    %[n_lo], %[t4]                \n\
    adcq    $0, %[t1]                     \n\
    movq    %[mmagic], %[n_lo]            \n\
    imulq   %[t6], %[n_lo]                \n\
    movq    %[n_lo], %%rdx                \n\
    mulxq   %[mod_lo], %[t0], %%rdx       \n\
    xorl    %k[n_hi], %k[n_hi]            \n\
    addq    %[t4], %[n_lo]                \n\
    setb    %b[n_hi]                      \n\
    addq    %[t6], %[t0]                  \n\
    adcq    %%rdx, %[n_lo]                \n\
    adcq    %[t1], %[n_hi]                \n\
":	"=&a"(n_lo), "=&r"(n_hi),[t0] "=&r"(t0),[t1] "=&r"(t1),[t2] "=&r"(t2),[t3] "=&r"(t3),[t4] "=&r"(t4),[t5] "=&r"(t5),
	    [t6] "=&r"(t6)
 :	    [n_lo] "0"(n_lo),[n_hi] "1"(n_hi),[mmagic] "r"(mmagic),[mod_lo] "r"(mod_lo)
 :	    "flags", "%rdx", "rcx");

	*res_lo = n_lo;
	*res_hi = n_hi;

#else
	ciosModSquare(res_lo, res_hi, mod_lo, mmagic);
	ciosModSquare(res_lo, res_hi, mod_lo, mmagic);
	ciosModSquare(res_lo, res_hi, mod_lo, mmagic);
#endif
}

bool optimizedSprpTest65(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);
#endif
	uint64_t bit;
	uint64_t res_lo, res_hi;
	uint64_t one_hi, one_lo;
	uint64_t scan_lo, k;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// count and trimm trailing zeroes
	k = my_ctz64(n_lo - 1);
	scan_lo = my_shr64(n_lo, k);

	// 1 in montgomery domain is 1 * 2^128 % m 
	// iterate over the 6 most significant bits
	// and enter montgomery demain
	uint64_t s = 5;
	s = (k <= 64 - s) ? s : 64 - k;
	uint64_t msbs = (s == 0) ? 1 : ((1 << s) + (n_lo >> (64 - s)));
	ciosConstants(msbs, n_lo, &one_lo, &one_hi, &res_lo, &res_hi);
	bit = (64 - s) - k;

	while (bit >= 3) {
		bit -= 3;
		// square and reduce 3 times
		ciosModSquare3(&res_lo, &res_hi, n_lo, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((scan_lo >> bit) & 7));
	}

	while (bit > 0) {
		bit -= 1;
		// square and reduce 3 times
		ciosModSquare(&res_lo, &res_hi, n_lo, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((scan_lo >> bit) & 1));
	}

	// make sure result is strictly less than the modulus
	ciosSubtract(&res_lo, &res_hi, n_lo);

	// if 1, return true
	if ((res_lo == one_lo) && (res_hi == one_hi))
		return true;

	// m-1 in montgomery domain is m - (1 * 2^128 % m)
	uint64_t m1_hi, m1_lo;
	uint8_t c;
	c = my_sbb64(0, n_lo, one_lo, &m1_lo);
	my_sbb64(c, 1, one_hi, &m1_hi);

	while (k > 1) {
		k -= 1;
		if ((res_lo == m1_lo) && (res_hi == m1_hi))
			return true;
		// square and reduce
		ciosModSquare(&res_lo, &res_hi, n_lo, mmagic);
		ciosSubtract(&res_lo, &res_hi, n_lo);
		if ((res_lo == one_lo) && (res_hi == one_hi))
			return false;

	}
	return ((res_lo == m1_lo) && (res_hi == m1_hi));
}

bool optimizedFermatTest65(uint64_t n_lo)
{
	uint64_t res_lo, res_hi, bit;
	uint64_t one_hi, one_lo;
	uint64_t scan_lo, k;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// 1 in montgomery domain is (1 * 2^128) % m
	// process the 6 most significant bits and convert to montgomery domain  
	uint64_t s = 5;
	uint64_t msbs = s == 0 ? 1 : ((1 << s) + (n_lo >> (64 - s)));
	ciosConstants(msbs, n_lo, &one_lo, &one_hi, &res_lo, &res_hi);
	bit = 64 - s;

	while (bit >= 5) {
		bit -= 3;
		// square and reduce 3 times
		ciosModSquare3(&res_lo, &res_hi, n_lo, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((n_lo >> bit) & 7));
	}

	while (bit > 1) {
		bit -= 1;
		// square and reduce
		ciosModSquare(&res_lo, &res_hi, n_lo, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((n_lo >> bit) & 1));
	}

	// make sure result is strictly less than the modulus
	ciosSubtract(&res_lo, &res_hi, n_lo);

	// - Euler's criterion 2^(n>>1) == legendre_symbol(2,n) (https://en.wikipedia.org/wiki/Euler%27s_criterion)
	// - Fermat primality check:
	//   (2^(n-1) == 1) mod n
	//
	// - Euler primality check:
	//   (2^(n>>1) == 1) mod n
	//   (2^(n>>1) == n-1) mod n
	uint64_t legendre = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	// when bit 1 and bit 2 are different, then n = 3 or 5 mod 8 and legendre(2,n) is -1, l = 1
	// when bit 1 and bit 2 are same, then n = 1 or 7 mod 8 and legendre(2,n) is 1, l = 0
	// compute n-1
	// check pseudo-primality
	//   return false : n is composite for sure
	//   return true : n is maybe prime, more tests are needed to prove primality 
	//
	uint64_t m1_hi, m1_lo;
	uint8_t c;
	// m-1 in montgomery domain is m - (1 * 2^128) % m
	c = my_sbb64(0, n_lo, one_lo, &m1_lo);
	my_sbb64(c, 1, one_hi, &m1_hi);

	return ((res_lo == (legendre ? m1_lo : one_lo)) && (res_hi == (legendre ? m1_hi : one_hi)));
}
