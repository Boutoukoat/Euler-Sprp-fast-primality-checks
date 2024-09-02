// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// fastest code
//
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>

#include "m64_utils.h"
#include "montgomery.h"
#include "optimized.h"
#include "m_reg.h"

typedef unsigned __int128 uint128_t;

// subtract the modulus 'mod' multiple times from the input number 'res', if needed
static inline void ciosSubtract128(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mod_hi)
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
	}
	while (b == 0);
	// get the saved values
	*res_lo = t_lo;
	*res_hi = t_hi;
}

// subtract the modulus i'mod' multiple times from the input number 'res', if needed
static inline void ciosSubtract64(uint64_t * res_lo, uint64_t mod_lo)
{
	uint64_t n_lo, n_hi;
	uint64_t t_lo, t_hi;
	uint8_t b;
	n_lo = *res_lo;
	// save, subtract the modulus until a borrows occurs
	do {
		t_lo = n_lo;
		b = my_sbb64(0, n_lo, mod_lo, &n_lo);
	}
	while (b == 0);
	// get the saved values
	*res_lo = t_lo;
}

static void ciosConstants64(uint64_t mod_lo, uint64_t * magic_lo)
{
	// computes 2^64 % mod
	uint64_t t = -mod_lo;	// 2^64-m
	t %= mod_lo;		// (2^64-m) % m
	*magic_lo = t;
#if PARANOID
	assert(*magic_lo <= mod_lo);
#endif
}

static void ciosConstants128(uint64_t mod_lo, uint64_t mod_hi, uint64_t * magic_lo, uint64_t * magic_hi)
{
	// computes 2^128 % mod
	uint128_t m = ((uint128_t) mod_hi << 64) + mod_lo;
	uint128_t t = -m;	// 2^128-m
	t %= m;			// (2^128-m) % m
	*magic_lo = (uint64_t) t;
	*magic_hi = (uint64_t) (t >> 64);
#if PARANOID
	assert(*magic_hi <= mod_hi);
#endif
}

static inline __attribute__((always_inline))
void ciosModSquare64(uint64_t * res_lo, uint64_t mod_lo, uint64_t mmagic)
{
	uint64_t n_lo = *res_lo;
	uint128_t cs, cc;
	uint64_t t0, t1, m, ignore;

	cc = (uint128_t) n_lo *n_lo;	// #1
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	t1 = (uint64_t) cc;
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t) m *mod_lo;	// #4
	cs += t0;
	cs = cs >> 64;
	cs += t1;
	t0 = (uint64_t) cs;
#if 0
	// not necessary with 2-bits guard
	cs = cs >> 64;
	t1 = (uint64_t) cs;
#endif
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	*res_lo = t0;
}

static inline __attribute__((always_inline))
void ciosModSquare128(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mod_hi, uint64_t mmagic)
{
	uint64_t n_lo = *res_lo, n_hi = *res_hi;
	uint128_t cs, cc;
	uint64_t t0, t1, t2, m, ignore;

	cc = (uint128_t) n_lo *n_lo;	// #1
	t0 = (uint64_t) cc;
	cc = cc >> 64;
	cc += (uint128_t) n_lo *(n_hi + n_hi);	// #2
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
	if (__builtin_constant_p(mod_hi) && mod_hi == 0) {
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 1) {
		cs += m;	//#3, removed at compile time
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 2) {
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 3) {
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
	} else {
		cs += (uint128_t) m *mod_hi;	//#3, removed at compile time
	}
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

	cc = (uint128_t) n_hi *n_hi;	// #6
	cc += t1;
	t1 = (uint64_t) cc;
	cc = cc >> 64;
	cc += t2;
	t2 = (uint64_t) cc;
#if 0
	// not necessary with 2-bits guard
	cc = cc >> 64;
	uint64_t t3 = (uint64_t) cc;
	assert(t3 == 0);
#endif
#if PARANOID
	assert(cc >> 64 == 0);
#endif

	m = t0 * mmagic;	// #3
	cs = (uint128_t) m *mod_lo;	// #8
	cs += t0;
	cs = cs >> 64;
	if (__builtin_constant_p(mod_hi) && mod_hi == 0) {
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 1) {
		cs += m;	//#3, removed at compile time
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 2) {
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
	} else if (__builtin_constant_p(mod_hi) && mod_hi == 3) {
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
		cs += m;	//#3, removed at compile time
	} else {
		cs += (uint128_t) m *mod_hi;	//#3, removed at compile time
	}
	cs += t1;
	t0 = (uint64_t) cs;
	cs = cs >> 64;
	cs += t2;
	t1 = (uint64_t) cs;
#if 0
	// not necessary with 2-bits guard
	cs = cs >> 64;
	cs += t3;
	t2 = (uint64_t) cs;
	assert(t2 == 0);
#endif
#if PARANOID
	assert(cs >> 64 == 0);
#endif

	*res_lo = t0;
	*res_hi = t1;

}

static inline __attribute__((always_inline))
void ciosModSquare3_128(uint64_t * res_lo, uint64_t * res_hi, uint64_t mod_lo, uint64_t mod_hi, uint64_t mmagic)
{
	ciosModSquare128(res_lo, res_hi, mod_lo, mod_hi, mmagic);
	ciosModSquare128(res_lo, res_hi, mod_lo, mod_hi, mmagic);
	ciosModSquare128(res_lo, res_hi, mod_lo, mod_hi, mmagic);
}

static inline __attribute__((always_inline))
void ciosModSquare3_64(uint64_t * res_lo, uint64_t mod_lo, uint64_t mmagic)
{
	ciosModSquare64(res_lo, mod_lo, mmagic);
	ciosModSquare64(res_lo, mod_lo, mmagic);
	ciosModSquare64(res_lo, mod_lo, mmagic);
}

static inline __attribute__((always_inline))
bool ciosFermatTest64(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);
#endif

	uint64_t res_lo, one_lo;
	int bit;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant 2^64 mod m
	ciosConstants64(n_lo, &one_lo);
	// constant (1 << msbits) * 2^64 mod m

	bit = 64 - my_clz64(n_lo);
	uint64_t msb_bits = bit < 5 ? bit - 1 : 3;
	uint64_t msb_mask = (1 << msb_bits) - 1;
	bit -= msb_bits;
	res_lo = one_lo << ((n_lo >> bit) & msb_mask);

	while (bit >= 5) {
		bit -= 3;
		// square and reduce
		ciosModSquare3_64(&res_lo, n_lo, mmagic);
		// shift
		res_lo <<= (n_lo >> bit) & 7;
	}

	while (bit > 1) {
		bit -= 1;
		// square and reduce
		ciosModSquare64(&res_lo, n_lo, mmagic);
		// shift
		res_lo <<= (n_lo >> bit) & 1;
	}

	// make sure result is strictly less than the modulus
	ciosSubtract64(&res_lo, n_lo);

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

	// check pseudo-primality
	//   return false : n is composite for sure
	//   return true : n is maybe prime, more tests are needed (like Lucas test)
	//
	return (res_lo == (legendre ? n_lo - one_lo : one_lo));
}

static inline __attribute__((always_inline))
bool ciosFermatTest128(uint64_t n_lo, uint64_t n_hi)
{
#if PARANOID
	assert((n_lo & 1) == 1);
#endif

	uint64_t res_lo, res_hi;
	uint64_t one_lo, one_hi;
	int bit;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant 2^128 mod m
	ciosConstants128(n_lo, n_hi, &one_lo, &one_hi);
	res_hi = one_hi;
	res_lo = one_lo;

	if (__builtin_constant_p(n_hi) && n_hi == 0) {
		bit = 64 - my_clz64(n_lo);
	uint64_t msb_bits = bit < 5 ? bit - 1 : 3;
	uint64_t msb_mask = (1 << msb_bits) - 1;
	bit -= msb_bits;
	my_shld64(&res_hi, &res_lo, (n_lo >> bit) & msb_mask);

	} else {
		if (n_hi == 0) {
			bit = 64 - my_clz64(n_lo);
	uint64_t msb_bits = bit < 5 ? bit - 1 : 3;
	uint64_t msb_mask = (1 << msb_bits) - 1;
	bit -= msb_bits;
	my_shld64(&res_hi, &res_lo, (n_lo >> bit) & msb_mask);

		} else {
			bit = 64 - my_clz64(n_hi);
			uint64_t msb_bits = bit < 4 ? bit : 3;
			uint64_t msb_mask = (1 << msb_bits) - 1;
			bit -= msb_bits;
			my_shld64(&res_hi, &res_lo, (n_hi >> bit) & msb_mask);

			while (bit >= 3) {
				bit -= 3;
				// square and reduce
				ciosModSquare3_128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
				// shift
				my_shld64(&res_hi, &res_lo, ((n_hi >> bit) & 7));
			}

			while (bit) {
				bit -= 1;
				// square and reduce
				ciosModSquare128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
				// shift
				my_shld64(&res_hi, &res_lo, ((n_hi >> bit) & 1));
			}

			bit = 64;
		}
	}
	while (bit >= 5) {
		bit -= 3;
		// square and reduce
		ciosModSquare3_128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((n_lo >> bit) & 7));
	}

	while (bit > 1) {
		bit -= 1;
		// square and reduce
		ciosModSquare128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((n_lo >> bit) & 1));
	}

	// make sure result is strictly less than the modulus
	ciosSubtract128(&res_lo, &res_hi, n_lo, n_hi);

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
	//
	// compute (m-1) * 2^128 % mod = modulus - one
	uint64_t m1_lo;
	uint64_t m1_hi;
	uint8_t c;
	c = my_sbb64(0, n_lo, one_lo, &m1_lo);
	my_sbb64(c, n_hi, one_hi, &m1_hi);
	// check pseudo-primality
	//   return false : n is composite for sure
	//   return true : n is maybe prime, more tests are needed (like Lucas test)
	//
	return ((res_lo == (legendre ? m1_lo : one_lo)) && (res_hi == (legendre ? m1_hi : one_hi)));
}

static inline __attribute__((always_inline))
bool ciosSprpTest64(uint64_t n_lo)
{
#if PARANOID
	assert((n_lo & 1) == 1);
#endif

	uint64_t res_lo;
	uint64_t one_lo;
	uint64_t s_lo;
	uint64_t m1_lo;
	uint64_t bit, k;
	uint8_t c;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant 1 * 2^64 mod m
	ciosConstants64(n_lo, &one_lo);

	// constant n-1 * 2^64
	m1_lo = n_lo - one_lo;

	// n - 1 = s * 2^k
	k = my_ctz64(n_lo - 1);
	s_lo = n_lo >> k;

	bit = 64 - my_clz64(s_lo);
	// constant (1 << msbits) * 2^64
	uint64_t msb_bits = bit < 5 ? bit - 1 : 3;
	uint64_t msb_mask = (1 << msb_bits) - 1;
	bit -= msb_bits;
	res_lo = one_lo << ((s_lo >> bit) & msb_mask);

	while (bit >= 5) {
		bit -= 3;
		// square and reduce
		ciosModSquare3_64(&res_lo, n_lo, mmagic);
		// shift
		res_lo <<= (s_lo >> bit) & 7;
	}

	while (bit > 0) {
		bit -= 1;
		// square and reduce
		ciosModSquare64(&res_lo, n_lo, mmagic);
		// shift
		res_lo <<= (s_lo >> bit) & 1;
	}

	ciosSubtract64(&res_lo, n_lo);
	// if res == 1 return true;
	if (res_lo == one_lo)
		return true;

	while (k > 1) {
		k -= 1;
		// if res == n-1 return true;
		if (res_lo == m1_lo)
			return true;
		// square and reduce
		ciosModSquare64(&res_lo, n_lo, mmagic);
		ciosSubtract64(&res_lo, n_lo);
		// if res == 1 return false;
		if (res_lo == one_lo)
			return false;
	}
	// if res == n-1 return true;
	return (res_lo == m1_lo);
}

static inline __attribute__((always_inline))
bool ciosSprpTest128(uint64_t n_lo, uint64_t n_hi)
{
#if PARANOID
	assert((n_lo & 1) == 1);
	assert((n_lo | n_hi) != 0);
#endif

	uint64_t res_lo, res_hi;
	uint64_t one_lo, one_hi;
	uint64_t s_lo, s_hi;
	uint64_t m1_lo, m1_hi;
	uint64_t bit, k;
	uint8_t c;
	// constant -1/m mod 2^64
	uint64_t mmagic = montgomeryInverse64(n_lo);

	// enter montgomery domain
	// constant one = 1 * 2^128 mod m
	ciosConstants128(n_lo, n_hi, &one_lo, &one_hi);

	// constant n-1 * 2^128 = n - one
	c = my_sbb64(0, n_lo, one_lo, &m1_lo);
	my_sbb64(c, n_hi, one_hi, &m1_hi);

	res_lo = one_lo;
	res_hi = one_hi;

	// n - 1 = s * 2^k
	k = my_ctz128(n_lo - 1, n_hi);
	s_lo = n_lo;
	s_hi = n_hi;
	my_shrd64(&s_hi, &s_lo, k);

	if (s_hi != 0) {
		bit = 64 - my_clz64(s_hi);
		uint64_t msb_bits = bit < 4 ? bit : 3;
		uint64_t msb_mask = (1 << msb_bits) - 1;
		bit -= msb_bits;
		my_shld64(&res_hi, &res_lo, (s_hi >> bit) & msb_mask);

		while (bit >= 3) {
			bit -= 3;
			// square and reduce
			ciosModSquare3_128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
			// shift
			my_shld64(&res_hi, &res_lo, ((s_hi >> bit) & 7));
		}

		while (bit > 0) {
			bit -= 1;
			// square and reduce
			ciosModSquare128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
			// shift
			my_shld64(&res_hi, &res_lo, ((s_hi >> bit) & 1));
		}

		bit = 64;
	} else {
		bit = 64 - my_clz64(s_lo);
		uint64_t msb_bits = bit < 5 ? bit - 1 : 3;
		uint64_t msb_mask = (1 << msb_bits) - 1;
		bit -= msb_bits;
		my_shld64(&res_hi, &res_lo, (s_lo >> bit) & msb_mask);
	}

	while (bit >= 3) {
		bit -= 3;
		// square and reduce
		ciosModSquare3_128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((s_lo >> bit) & 7));
	}

	while (bit > 0) {
		bit -= 1;
		// square and reduce
		ciosModSquare128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
		// shift
		my_shld64(&res_hi, &res_lo, ((s_lo >> bit) & 1));
	}

	// Make sure the result is less than the modulus
	ciosSubtract128(&res_lo, &res_hi, n_lo, n_hi);

	// if res == 1 return true;
	if ((res_lo == one_lo) && (res_hi == one_hi))
		return true;

	while (k > 1) {
		k -= 1;
		// if res == n-1 return true;
		if ((res_lo == m1_lo) && (res_hi == m1_hi))
			return true;
		// square and reduce
		ciosModSquare128(&res_lo, &res_hi, n_lo, n_hi, mmagic);
		ciosSubtract128(&res_lo, &res_hi, n_lo, n_hi);
	// if res == 1 return false;
	if ((res_lo == one_lo) && (res_hi == one_hi))
		return false;

	}
	// if res == n-1 return true;
	return ((res_lo == m1_lo) && (res_hi == m1_hi));
}

// ------------------------------------------------------------
// Fermat + Euler primality test base 2
// ------------------------------------------------------------

bool optimizedFermatTest(uint64_t n_lo, uint64_t n_hi)
{
	// dispatch to the fastest code for each bit size
	switch (n_hi) {
	case 1:
		return optimizedFermatTest65(n_lo);	// 65 bits
	case 2:
		return ciosFermatTest128(n_lo, 2);	// 66 bits 
	case 3:
		return ciosFermatTest128(n_lo, 3);	// 66 bits
	case 0:
		if (n_lo >> 8 == 0) {
			return montgomeryFermatTest8(n_lo);	// from 0 to 8 bits (blitz version)
		}
		if (n_lo >> 55 != 0) {
			return ciosFermatTest128(n_lo, 0);	// 63-64 bits
		}
		return ciosFermatTest64(n_lo);	// from 3 to 62 bits
	default:
		if (n_hi >> 55) {
			// there are shifts without reduction up to 7 bits in the optimized code
			return montgomeryFermatTest128(n_lo, n_hi);	// 120-127 bits
		}
		return ciosFermatTest128(n_lo, n_hi);	// 67-119 bits
	}
}

// ------------------------------------------------------------
// SPRP primality test base 2
// ------------------------------------------------------------

bool optimizedSprpTest(uint64_t n_lo, uint64_t n_hi)
{
	// dispatch to the fastest code for each bit size
	switch (n_hi) {
	case 1:
		return optimizedSprpTest65(n_lo);	// 65 bits
	case 2:
		return ciosSprpTest128(n_lo, 2);	// 66 bits 
	case 3:
		return ciosSprpTest128(n_lo, 3);	// 66 bits
	case 0:
		if (n_lo >> 8 == 0) {
			return montgomeryFermatTest8(n_lo);	// from 0 to 8 bits (blitz version)
		}
		if (n_lo >> 55 != 0) {
			// there are shifts without reduction up to 7 bits in the optimized code.
			// 64 bits - 7 bit shift - 2 guard bits = 55 bits.
			return ciosSprpTest128(n_lo, 0);	// from 56 to 64 bits.
		}
		return ciosSprpTest64(n_lo);	// 3-62 bits
	default:
		if (n_hi >> 55) {
			// there are shifts without reduction up to 7 bits in the optimized code.
			// 128 bits - 7 bit shift - 2 guard bits = 119 bits.
			// this requires code where shifts do not overflow
			return montgomerySprpTest128(n_lo, n_hi);	// 120-127 bits
		}
		return ciosSprpTest128(n_lo, n_hi);	// 67-119 bits
	}
}
