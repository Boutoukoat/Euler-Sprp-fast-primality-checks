// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// moderately optimized primality tests
//
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "m64_utils.h"
#include "montgomery.h"
#include "m_reg.h"

typedef unsigned __int128 uint128_t;

// subtract until the result is less than the modulus
static void ciosSubtract64(uint64_t * res_lo, uint64_t carries, uint64_t mod_lo)
{
	uint64_t n_lo;
	uint64_t t_lo;
	uint8_t b;
	n_lo = *res_lo;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
		if (__builtin_constant_p(carries) && carries == 0) {
		} else {
			b = my_sbb64(b, carries, 0, &carries);
		}
	}
	while (b == 0);
	// get the saved value when a borrow occurs
	*res_lo = t_lo;
}

// subtract until the result is less than the modulus
static void ciosSubtract128(uint64_t * res_lo, uint64_t * res_hi, uint64_t carries, uint64_t mod_lo, uint64_t mod_hi)
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
		b = my_sbb64(b, n_hi, mod_hi, &n_hi);
		if (__builtin_constant_p(carries) && carries == 0) {
		} else {
			b = my_sbb64(b, carries, 0, &carries);
		}
	}
	while (b == 0);
	// get the saved values when a borrow occurs
	*res_lo = t_lo;
	*res_hi = t_hi;
}

// shift and subtract until the bits in eccess fade out, result is possibly > modulus
static void ciosModShift64(uint64_t * res_lo, uint64_t mod_lo, uint64_t shift)
{
#if PARANOID
	//  this code assume that the shift amount is only a few bits
	assert(shift < 64);
#endif
	uint64_t t0 = *res_lo, t1 = 0;
	uint8_t c;

	my_shld64(&t1, &t0, shift);

	while (t1) {
		c = my_sbb64(0, t0, mod_lo, &t0);
		my_sbb64(c, t1, 0, &t1);
	}
	*res_lo = t0;
}

// shift and subtract until the bits in eccess fade out, result is possibly > modulus
static void ciosModShift128(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mod_hi, uint64_t shift)
{
#if PARANOID
	//  this code assume that the shift amount is only a few bits
	assert(shift < 64);
#endif
	uint64_t t0 = *res_lo, t1 = *res_hi, t2 = 0;
	uint8_t c;

	my_shld64(&t2, &t1, &t0, shift);

	while (t2) {
		c = my_sbb64(0, t0, mod_lo, &t0);
		c = my_sbb64(c, t1, mod_hi, &t1);
		my_sbb64(c, t2, 0, &t2);
	}
	*res_lo = t0;
	*res_hi = t1;
}

// computes 2^64 %  mod
static void ciosConstants64(uint64_t mod_lo, uint64_t * magic2_lo)
{
	// computes 2^64 %  mod
	uint128_t t = (uint128_t) 1 << 64;
	t %= mod_lo;
	*magic2_lo = (uint64_t) t;
}

// computes 2^128 %  mod
static void ciosConstants128(uint64_t mod_lo, uint64_t mod_hi, uint64_t * magic2_lo, uint64_t * magic2_hi)
{
	// computes 2^128 % mod = 2 * (2^127) % mod
	uint128_t mod = ((uint128_t) mod_hi << 64) + mod_lo;

	// step 1  1 << 127 % mod
	uint128_t t = (uint128_t) 1 << 127;
	t %= mod;

	// step 2  ((1 << 127 % mod) << 1) % mod
	uint64_t t0 = (uint64_t) t;
	uint64_t t1 = (uint64_t) (t >> 64);
	uint64_t t2 = 0;
	ciosModShift128(&t0, &t1, mod_lo, mod_hi, 1);
	*magic2_lo = t0;
	*magic2_hi = t1;
}

// Montgomery modular multiplication    res = (res * b) % mod
static void ciosModMul64(uint64_t * res_lo, uint64_t b_lo, uint64_t mod_lo, uint64_t mmagic)
{
	uint64_t a_lo = *res_lo;
	uint128_t cs, cc;
	uint64_t t0, t1, m;

	cc = (uint128_t) a_lo *b_lo;	// #1
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	t1 = (uint64_t) cc;

	m = t0 * mmagic;	// #2
	cs = (uint128_t) m *mod_lo;	// #3
	cs += t0;
	cs = cs >> 64;
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	t1 = (uint64_t) cs;
	if (t1) {
		ciosSubtract64(&t0, t1, mod_lo);
	}

	*res_lo = t0;
}

static void ciosModMul128(uint64_t * res_lo, uint64_t * res_hi, uint64_t b_lo, uint64_t b_hi, uint64_t mod_lo, uint64_t mod_hi,
			  uint64_t mmagic)
{
	uint64_t a_lo = *res_lo, a_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, t3, m, ignore;

	cc = (uint128_t) a_lo *b_lo;	// #1
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	cc += (uint128_t) a_lo *b_hi;	// #2
	t1 = (uint64_t) cc;
	cc = cc >> 64;
	t2 = (uint64_t) cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t) m *mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t) m *mod_hi;	// #5
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t) cs;
	cs = cs >> 64;
	t2 = (uint64_t) cs;
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	cc = (uint128_t) a_hi *b_lo;	// #6
	cc += t0;
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	cc += (uint128_t) a_hi *b_hi;	// #7
	cc += t1;
	t1 = (uint64_t) cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t) cc;
	cc = cc >> 64;
	t3 = (uint64_t) cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #8
	cs = (uint128_t) m *mod_lo;	// #9
	cs += t0;
	cs = cs >> 64;
	cs += (uint128_t) m *mod_hi;	// #10
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t) cs;
	if (t2) {
		ciosSubtract128(&t0, &t1, t2, mod_lo, mod_hi);
	}

	*res_lo = t0;
	*res_hi = t1;
}

bool montgomeryFermatTest8(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);	// odd number <= 255
	assert(n_lo >> 8 == 0);	// for very small numbers
#endif
	// 2^((n_lo-1)/2) mod n_lo
	uint64_t res64;
	if (n_lo <= 127) {
		res64 = 1;
		res64 <<= (n_lo >> 1);	// (n_lo >> 1) is <= 63
		res64 %= n_lo;
	} else {
		uint128_t res128 = 1;
		res128 <<= (n_lo >> 1);	// (n_lo >> 1) is <= 127
		res128 %= n_lo;
		res64 = (uint64_t) res128;
	}
	// Euler criterion
	uint64_t l = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	return (res64 == (l ? n_lo - 1 : 1));
}

bool montgomeryFermatTest64(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);	// odd number
#endif
	uint64_t res_lo;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant 2^64 mod m
	ciosConstants64(n_lo, &res_lo);
	// value 2 * 2^64
	ciosModShift64(&res_lo, n_lo, 1);

	int bit = 63 - my_clz64(n_lo);
	while (bit >= 3) {
		bit -= 2;
		// square and reduce
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		// shift and reduce
		ciosModShift64(&res_lo, n_lo, ((n_lo >> bit) & 3));
	}
	while (bit > 1) {
		bit -= 1;
		// square and reduce
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		// shift and reduce
		ciosModShift64(&res_lo, n_lo, ((n_lo >> bit) & 1));
	}
	// exit montgomery domain
	ciosModMul64(&res_lo, 1, n_lo, mmagic);
	// make sure result is < modulus
	ciosSubtract64(&res_lo, 0, n_lo);

	// - Euler's criterion 2^(n>>1) == legendre_symbol(2,n) (https://en.wikipedia.org/wiki/Euler%27s_criterion)
	//
	// - Euler primality check:
	//   (2^(n>>1) == 1) mod n
	//   (2^(n>>1) == n-1) mod n
	uint64_t l = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	// when bit 1 and bit 2 are different, then n = 3 or 5 mod 8 and legendre(2,n) is -1, l = 1
	// when bit 1 and bit 2 are same, then n = 1 or 7 mod 8 and legendre(2,n) is 1, l = 0

	// check pseudo-primality
	//   return false : n is composite for sure
	//   return true : n is maybe prime, more tests are needed (like Lucas test)
	//
	return (res_lo == (l ? n_lo - 1 : 1));
}

bool montgomeryFermatTest128(uint64_t n_lo, uint64_t n_hi)
{
#if PARANOID
	assert((n_lo & 1) == 1);	// odd number
	assert(n_hi > 0);	// large number
#endif

	uint64_t res_lo, res_hi;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant 2^128 mod m
	ciosConstants128(n_lo, n_hi, &res_lo, &res_hi);
	// constant 2 * 2^128
	ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, 1);

	int bit = 63 - my_clz64(n_hi);
	while (bit >= 2) {
		bit -= 2;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift and reduce
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((n_hi >> bit) & 3));
	}
	while (bit) {
		bit -= 1;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift and reduce
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((n_hi >> bit) & 1));
	}

	bit = 64;
	while (bit > 2) {
		bit -= 2;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((n_lo >> bit) & 3));
	}
	while (bit > 1) {
		bit -= 1;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((n_lo >> bit) & 1));
	}

	// exit montgomery domain
	ciosModMul128(&res_lo, &res_hi, 1, 0, n_lo, n_hi, mmagic);
	// make sure result is < modulus
	ciosSubtract128(&res_lo, &res_hi, 0, n_lo, n_hi);

	// - Euler's criterion 2^(n>>1) == legendre_symbol(2,n) (https://en.wikipedia.org/wiki/Euler%27s_criterion)
	//
	// - Euler primality check:
	//   (2^(n>>1) == 1) mod n
	//   (2^(n>>1) == n-1) mod n
	uint64_t l = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	// when bit 1 and bit 2 are different, then n = 3 or 5 mod 8 and legendre(2,n) is -1, l = 1
	// when bit 1 and bit 2 are same, then n = 1 or 7 mod 8 and legendre(2,n) is 1, l = 0
	// check pseudo-primality
	//   return false : n is composite for sure
	//   return true : n is maybe prime, more tests are needed (like Lucas test)
	//
	return ((res_lo == (l ? n_lo - 1 : 1)) && (res_hi == (l ? n_hi : 0)));
}

bool montgomerySprpTest64(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);	// odd number
	assert(n_lo > 3);	// troubles in tests < 4
#endif
	uint64_t res_lo, one_lo, s_lo, m1_lo;
	uint64_t bit, k;
	uint8_t c;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// constant 1 * 2^64 mod m
	ciosConstants64(n_lo, &one_lo);

	// constant n-1 * 2^64
	m1_lo = n_lo - one_lo;

	// constant 2 * 2^64
	res_lo = one_lo;
	ciosModShift64(&res_lo, n_lo, 1);

	// n - 1 = s * 2^k
	k = my_ctz64(n_lo - 1);
	s_lo = my_shr64(n_lo, k);

	bit = 63 - my_clz64(s_lo);

	while (bit > 2) {
		bit -= 2;
		// square and reduce
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		// shift
		ciosModShift64(&res_lo, n_lo, (s_lo >> bit) & 3);
	}

	while (bit > 0) {
		bit -= 1;
		// square and reduce
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		// shift
		ciosModShift64(&res_lo, n_lo, (s_lo >> bit) & 1);
	}

	ciosSubtract64(&res_lo, 0, n_lo);
	// if res == 1 return true;
	if (res_lo == one_lo)
		return true;

	while (k > 1) {
		k -= 1;
		// if res == n-1 return true;
		if (res_lo == m1_lo)
			return true;
		// square and reduce
		ciosModMul64(&res_lo, res_lo, n_lo, mmagic);
		ciosSubtract64(&res_lo, 0, n_lo);
		if (res_lo == one_lo)
			return false;
	}
	// if res == n-1 return true;
	return (res_lo == m1_lo);
}

bool montgomerySprpTest128(uint64_t n_lo, uint64_t n_hi)
{
#if PARANOID
	assert((n_lo & 1) == 1);	// odd number
#endif
	uint64_t res_lo, res_hi;
	uint64_t one_lo, one_hi;
	uint64_t s_lo, s_hi;
	uint64_t m1_lo, m1_hi;
	uint64_t bit, k;
	uint8_t c;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// constant 1 * 2^128 mod m
	ciosConstants128(n_lo, n_hi, &one_lo, &one_hi);
	ciosSubtract128(&one_lo, &one_hi, 0, n_lo, n_hi);

	// constant n-1 * 2^128
	c = my_sbb64(0, n_lo, one_lo, &m1_lo);
	my_sbb64(c, n_hi, one_hi, &m1_hi);

	// constant 2 * 2^128
	res_lo = one_lo;
	res_hi = one_hi;
	ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, 1);

	// n - 1 = s * 2^k
	k = my_ctz128(n_lo - 1, n_hi);
	s_lo = n_lo;
	s_hi = n_hi;
	my_shrd64(&s_hi, &s_lo, k);

	if (s_hi != 0) {
		bit = 63 - my_clz64(s_hi);
		while (bit > 2) {
			bit -= 2;
			// square and reduce
			ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
			ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
			// shift
			ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((s_hi >> bit) & 3));
		}

		while (bit > 0) {
			bit -= 1;
			// square and reduce
			ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
			// shift
			ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((s_hi >> bit) & 1));
		}

		bit = 64;
	} else {
		bit = 63 - my_clz64(s_lo);
	}

	while (bit > 2) {
		bit -= 2;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((s_lo >> bit) & 3));
	}

	while (bit > 0) {
		bit -= 1;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		// shift
		ciosModShift128(&res_lo, &res_hi, n_lo, n_hi, ((s_lo >> bit) & 1));
	}

	// if res == 1 return true;
	ciosSubtract128(&res_lo, &res_hi, 0, n_lo, n_hi);
	if ((res_lo == one_lo) && (res_hi == one_hi))
		return true;

	while (k > 1) {
		k -= 1;
		// if res == n-1 return true;
		if ((res_lo == m1_lo) && (res_hi == m1_hi))
			return true;
		// square and reduce
		ciosModMul128(&res_lo, &res_hi, res_lo, res_hi, n_lo, n_hi, mmagic);
		ciosSubtract128(&res_lo, &res_hi, 0, n_lo, n_hi);
		if ((res_lo == one_lo) && (res_hi == one_hi))
			return false;

	}
	// if res == n-1 return true;
	return ((res_lo == m1_lo) && (res_hi == m1_hi));
}

bool montgomeryFermatTest(uint64_t n_lo, uint64_t n_hi)
{
	// dispatch to textbook code
	if (n_hi == 0) {
		if (n_lo >> 2 == 0) {
			return montgomeryFermatTest8(n_lo);	// handle 2 bit modulus
		}
		return montgomeryFermatTest64(n_lo);	// <= 64 bits  (textbook version)
	}
	return montgomeryFermatTest128(n_lo, n_hi);	// > 64 bits  (textbook version)
}

bool montgomerySprpTest(uint64_t n_lo, uint64_t n_hi)
{
	// dispatch to textbook code
	if (n_hi == 0) {
		if (n_lo >> 2 == 0) {
			return montgomeryFermatTest8(n_lo);	// handle 2 bit modulus
		}
		return montgomerySprpTest64(n_lo);	// <= 64 bits  (textbook version)
	}
	return montgomerySprpTest128(n_lo, n_hi);	// > 64 bits  (textbook version)
}
