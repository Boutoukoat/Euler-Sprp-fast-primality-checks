#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "m128_utils.h"
#include "barrett.h"
#include "m_reg.h"

// t2::t1::t0 -= s2::s1::s0
static uint64_t barrett_sub64(uint64_t * t0, uint64_t * t1, uint64_t s0, uint64_t s1)
{
	uint64_t r0, r1, tmp, borrow;
	r0 = *t0 - s0;
	borrow = r0 > *t0;
	tmp = *t1 - borrow;
	borrow = tmp > *t1;
	r1 = tmp - s1;
	borrow |= r1 > tmp;
	*t0 = r0;
	*t1 = r1;
	return borrow;
}

// t2::t1::t0 -= s2::s1::s0
static uint64_t barrett_sub(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2, tmp, borrow;
	r0 = *t0 - s0;
	borrow = r0 > *t0;
	tmp = *t1 - borrow;
	borrow = tmp > *t1;
	r1 = tmp - s1;
	borrow |= r1 > tmp;
	tmp = *t2 - borrow;
	borrow = tmp > *t2;
	r2 = tmp - s2;
	borrow |= r2 > tmp;
	*t0 = r0;
	*t1 = r1;
	*t2 = r2;
	return borrow;
}

// t2::t1::t0 += s2::s1::s0
static uint64_t barrett_add(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2, tmp, carry;
	r0 = *t0 + s0;
	carry = r0 < s0;
	tmp = *t1 + carry;
	carry = tmp < carry;
	r1 = tmp + s1;
	carry |= r1 < s1;
	tmp = *t2 + carry;
	carry = tmp < carry;
	r2 = tmp + s2;
	carry |= r2 < s2;
	*t0 = r0;
	*t1 = r1;
	*t2 = r2;
	return carry;
}

// t2::t1::t0 >>= s
static void barrett_shr(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s)
{
	uint64_t r0, r1, r2;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	while (s >= 64) {
		r0 = r1;
		r1 = r2;
		r2 = 0;
		s -= 64;
	}
	if (s) {
		r0 = (r0 >> s) | (r1 << (64 - s));
		r1 = (r1 >> s) | (r2 << (64 - s));
		r2 = (r2 >> s);
	}
	*t0 = r0;
	*t1 = r1;
	*t2 = r2;
}

// t2::t1::t0 <<= s
static void barrett_shl(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s)
{
	uint64_t r0, r1, r2;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	while (s >= 64) {
		r2 = r1;
		r1 = r0;
		r0 = 0;
		s -= 64;
	}
	if (s) {
		r2 = (r2 << s) | (r1 >> (64 - s));
		r1 = (r1 << s) | (r0 >> (64 - s));
		r0 = (r0 << s);
	}
	*t0 = r0;
	*t1 = r1;
	*t2 = r2;
}

// t2::t1::t0 < s2::s1:s0 ?
static uint64_t barrett_cmp(uint64_t t0, uint64_t t1, uint64_t t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2, tmp, borrow;
	r0 = t0 - s0;
	borrow = r0 > t0;
	tmp = t1 - borrow;
	borrow = tmp > t1;
	r1 = tmp - s1;
	borrow |= r1 > tmp;
	tmp = t2 - borrow;
	borrow = tmp > t2;
	r2 = tmp - s2;
	borrow |= r2 > tmp;
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

// conditional subtract
// t2::t1::t0 < s2::s1:s0 ? t : t - s
static uint64_t barrett_cmp_sub(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2, borrow;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	borrow = barrett_sub(&r0, &r1, &r2, s0, s1, s2);
	if (!borrow) {
		*t0 = r0;
		*t1 = r1;
		*t2 = r2;
	}
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

// conditional subtract
// t1::t0 < s1:s0 ? t : t - s
static uint64_t barrett_cmp_sub64(uint64_t * t0, uint64_t * t1, uint64_t s0, uint64_t s1)
{
	uint64_t r0, r1, borrow;
	r0 = *t0;
	r1 = *t1;
	borrow = barrett_sub64(&r0, &r1, s0, s1);
	if (!borrow) {
		*t0 = r0;
		*t1 = r1;
	}
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

static void barrett_reduce64(uint64_t * t0, uint64_t * t1, uint64_t barrett_lo, uint64_t barrett_hi, uint64_t mod_lo)
{
	// process Barrett coefficient
	uint128_t ee_hi, ee_lo, cc, ss0, ss1;
	uint64_t e_hi, e_lo, s0, s1, borrow;
	s0 = *t0;
	s1 = *t1;
	ee_hi = (uint128_t) barrett_hi *s1;
	ee_lo = (uint128_t) barrett_hi *s0;
	ee_hi += (ee_lo >> 64);
	ee_lo = (uint128_t) barrett_lo *s1;
	ee_hi += (ee_lo >> 64);
	// keep the most significant bits of Barrett estimate
	e_hi = (uint64_t) (ee_hi >> 64);
	e_lo = (uint64_t) (ee_hi);

	// multiply by the modulus
	cc = (uint128_t) e_lo *mod_lo;
	s0 = (uint64_t) cc;
	cc >>= 64;
	cc += (uint128_t) e_hi *mod_lo;
	s1 = (uint64_t) cc;

#if PARANOID
	assert(cc >> 64 == 0);
#endif

	// subtract from initial number, by Barrett construction cannot underflow.
	borrow = barrett_sub64(t0, t1, s0, s1);
#if PARANOID
	assert(borrow == 0);
#endif
	while (*t1) {
		barrett_sub64(t0, t1, mod_lo, 0);
	}
}

static void barrett_reduce(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t barrett_lo, uint64_t barrett_hi, uint64_t mod_lo,
			   uint64_t mod_hi)
{
	// process Barrett coefficient
	uint128_t ee_hi, ee_lo, cc, ss0, ss1, ss2;
	uint64_t e_hi, e_lo, s0, s1, s2, borrow;
	s1 = *t1;
	s2 = *t2;
	ee_hi = (uint128_t) barrett_hi *s2;
	ee_lo = (uint128_t) barrett_hi *s1;
	ee_hi += (ee_lo >> 64);
	ee_lo = (uint128_t) barrett_lo *s2;
	ee_hi += (ee_lo >> 64);
	// keep the most significant bits of Barrett estimate
	e_hi = (uint64_t) (ee_hi >> 64);
	e_lo = (uint64_t) (ee_hi);

	// multiply by the modulus
	ss2 = (uint128_t) e_hi *mod_hi;
	ss1 = (uint128_t) e_hi *mod_lo;
	ss2 += (ss1 >> 64);
	ss1 = (uint64_t) ss1;
	cc = (uint128_t) e_lo *mod_hi;
	ss2 += (cc >> 64);
	ss1 += (uint64_t) cc;
	ss0 = (uint128_t) e_lo *mod_lo;

	cc = ss0;
	s0 = (uint64_t) cc;
	cc >>= 64;
	cc += ss1;
	s1 = (uint64_t) cc;
	cc >>= 64;
	cc += ss2;
	s2 = (uint64_t) cc;

	// subtract from initial number, by Barrett construction cannot underflow.
	borrow = barrett_sub(t0, t1, t2, s0, s1, s2);
#if PARANOID
	assert(borrow == 0);
#endif
	while (*t2) {
		barrett_sub(t0, t1, t2, mod_lo, mod_hi, 0);
	}
}

static void barrett_square(uint128_t v, uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t * t3)
{
	uint128_t tt0, tt1, tt2, tt3;
	uint64_t v0, v1;

	// t = v ^ 2
	v0 = (uint64_t) v;
	v1 = (uint64_t) (v >> 64);
	tt0 = (uint128_t) v0 *v0;
	tt1 = (uint128_t) v1 *v0;
	tt2 = (uint128_t) v1 *v1;

	// propagate carries along multiplication steps without overflow
	tt3 = tt2 >> 64;
	tt2 = (uint64_t) tt2;
	tt2 += tt1 >> 64;
	tt2 += tt1 >> 64;
	tt1 = (uint64_t) tt1;
	tt1 += tt1;
	tt1 += tt0 >> 64;
	tt0 = (uint64_t) tt0;

	// make t0, t1, t2, t3 64 bit numbers;
	*t0 = (uint64_t) tt0;
	tt2 += tt1 >> 64;
	*t1 = (uint64_t) tt1;
	tt3 += tt2 >> 64;
	*t2 = (uint64_t) tt2;
	*t3 = (uint64_t) tt3;
}

// compute 2^192 / modulus
static void reciprocal(uint64_t mod_lo, uint64_t mod_hi, uint64_t * barrett_lo, uint64_t * barrett_hi)
{
#if PARANOID
	assert(mod_hi > 0);
#endif

	// Variables used to correct the estimated quotient.
	uint128_t q, r, c1, c2;
	bool c3, c4;

	uint64_t shift = my_clz64(mod_hi);
	my_shld64(&mod_hi, &mod_lo, shift);
	uint128_t mod = ((uint128_t) mod_hi << 64) + mod_lo;
	uint128_t numhi = (uint128_t) 1 << (64 + shift);

	// We wish to compute q1 = [n3 n2 n1] / [d1 d0].
	// Estimate q1 as [n3 n2] / [d1], and then correct it.
	// Note while q may be 2 digits, q1 is always 1 digit.
	q = numhi / mod_hi;
	r = numhi - q * mod_hi;
	c1 = q * mod_lo;
	c2 = r << 64;
	c3 = c1 > c2;
	c4 = (c1 - c2) > mod;
	q -= c3 ? (c4 ? 2 : 1) : 0;

	// Compute the true (partial) remainder.
	uint128_t rem = (numhi << 64) - mod * (uint64_t) q;
	*barrett_hi = (uint64_t) q;

	// We wish to compute q0 = [rem1 rem0 n0] / [d1 d0].
	// Estimate q0 as [rem1 rem0] / [d1] and correct it.
	q = rem / mod_hi;
	r = rem - q * mod_hi;
	c1 = q * mod_lo;
	c2 = r << 64;
	c3 = c1 > c2;
	c4 = (c1 - c2) > mod;
	q -= c3 ? (c4 ? 2 : 1) : 0;
	*barrett_lo = (uint64_t) q;
}

// compute 2^128 / modulus
static void reciprocal64(uint64_t mod_lo, uint64_t * barrett_lo, uint64_t * barrett_hi)
{
#if PARANOID
	assert(mod_lo > 0);
#endif

	uint128_t q, r, c1, c2;
	bool c3, c4;

	uint64_t shift = 64 + my_clz64(mod_lo);
	uint128_t mod = (uint128_t) mod_lo << shift;
	uint64_t mod_hi = (uint64_t) (mod >> 64);
	uint128_t numhi = (uint128_t) 1 << shift;

	// We wish to compute q1 = [n3 n2 n1] / [d1 d0].
	// Estimate q1 as [n3 n2] / [d1], and then correct it.
	// Note while qhat may be 2 digits, q1 is always 1 digit.
	q = numhi / mod_hi;
	r = numhi - q * mod_hi;
	c1 = 0;
	c2 = r << 64;
	c3 = c1 > c2;
	c4 = (c1 - c2) > mod;
	q -= c3 ? (c4 ? 2 : 1) : 0;
	uint128_t rem = (numhi << 64) - mod * (uint64_t) q;
	*barrett_hi = (uint64_t) q;

	// We wish to compute q0 = [rem1 rem0 n0] / [d1 d0].
	// Estimate q0 as [rem1 rem0] / [d1] and correct it.
	q = rem / mod_hi;
	r = rem - q * mod_hi;
	c1 = 0;
	c2 = r << 64;
	c3 = c1 > c2;
	c4 = (c1 - c2) > mod;
	q -= c3 ? (c4 ? 2 : 1) : 0;
	*barrett_lo = (uint64_t) q;
}

static uint128_t barrett_mod_square(uint128_t v, uint64_t mod_lo, uint64_t mod_hi, uint64_t barrett_lo, uint64_t barrett_hi)
{
#if PARANOID
	assert(mod_hi > 0);	// at least bit modulus
	assert(mod_lo & 1);	// odd modulus
#endif
	// square first, reduce after
	uint64_t t0, t1, t2, t3, borrow;
	barrett_square(v, &t0, &t1, &t2, &t3);

	// modular reduction
	// printf("r2  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 
	barrett_reduce(&t1, &t2, &t3, barrett_lo, barrett_hi, mod_lo, mod_hi);
	// printf("r1  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 
	barrett_reduce(&t0, &t1, &t2, barrett_lo, barrett_hi, mod_lo, mod_hi);
	// printf("r0  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 

	// subtract the modulus until the result is less than modulus
	// this loop should run 0, 1 or 2 times, not much more (t < 2 * modulus for sure)
	do {
		// subtract. when a borrow occurs, it was too late.
		borrow = barrett_cmp_sub(&t0, &t1, &t2, mod_lo, mod_hi, 0);
	}
	while (!borrow);

	v = ((uint128_t) t1 << 64) + t0;
	return v;
}

static uint64_t barrett_mod_square64(uint64_t v, uint64_t mod_lo, uint64_t barrett_lo, uint64_t barrett_hi)
{
#if PARANOID
	assert(mod_lo > 0);	// non_zero modulus
#endif
	// square first, reduce after
	uint64_t t0, t1, borrow;
	uint128_t t = (uint128_t) v * v;
	t0 = (uint64_t) t;
	t1 = (uint64_t) (t >> 64);

	// modular reduction
	// printf("r2  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 
	barrett_reduce64(&t0, &t1, barrett_lo, barrett_hi, mod_lo);
	// printf("r0  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 

	// subtract the modulus until the result is less than modulus
	// this loop should run 0, 1 or 2 times, not much more (t < 2 * modulus for sure)
	do {
		// subtract. when a borrow occurs, it was too late.
		borrow = barrett_cmp_sub64(&t0, &t1, mod_lo, 0);
	}
	while (!borrow);

	return t0;
}

// Fermat test with Euler criterion.
static bool barrettFermatTest128(uint64_t n_lo, uint64_t n_hi)
{
	uint128_t n = ((uint128_t) n_hi << 64) + n_lo;

#if PARANOID
	// assume n odd
	assert((n & 1) == 1);
#endif
	uint64_t barrett_lo, barrett_hi;
	reciprocal(n_lo, n_hi, &barrett_lo, &barrett_hi);

	uint128_t result = 1;
	int bit = 128;
	// find  the msb
	while (bit > 1) {
		bit--;
		if ((n >> bit) & 1) {
			result = 2;
			break;
		}
	}
	while (bit > 1) {
		bit--;
		// square
		result = barrett_mod_square(result, n_lo, n_hi, barrett_lo, barrett_hi);

		if ((n >> bit) & 1) {
			// multiply by 2
			uint128_t t = result;
			result += result;
			if (result < t) {
				// wrap-around occured , this could occur if n is 127 bits.
				// input condition : result < n  then  2*result < 2*n
				t = result;
				result = t - n;	// subtract the modulus, ignore carry and borrow.
#if PARANOID
				// already checked carry was 1 and now check borrow is 1
				assert(result >= t);	// another wrap-around must occur
				assert(result < n);	// 2*result - n < n
#endif
			}
			while (result >= n) {
				result -= n;
			}
		}
	}
	// Euler criterion
	uint64_t legendre = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	return result == (legendre ? n - 1 : 1);
}

// Fermat test with Euler criterion.
static bool barrettFermatTest64(uint64_t n_lo)
{
#if PARANOID
	// assume n odd
	assert((n_lo & 1) == 1);
#endif
	uint64_t barrett_lo, barrett_hi;
	reciprocal64(n_lo, &barrett_lo, &barrett_hi);

	uint64_t result = 1;
	int bit = 64;
	// find  the msb
	while (bit > 1) {
		bit--;
		if ((n_lo >> bit) & 1) {
			result = 2;
			break;
		}
	}
	while (bit > 1) {
		bit--;
		// square
		result = barrett_mod_square64(result, n_lo, barrett_lo, barrett_hi);

		if ((n_lo >> bit) & 1) {
			// multiply by 2
			uint128_t t = result;
			t += result;
			while (t >= n_lo) {
				t -= n_lo;
			}
			result = (uint64_t) t;
		}
	}
	// Euler criterion
	uint64_t legendre = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	return result == (legendre ? n_lo - 1 : 1);
}

// Fermat test with Euler criterion.
bool barrettFermatTest(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return barrettFermatTest128(n_lo, n_hi);
	}
	return barrettFermatTest64(n_lo);
}
