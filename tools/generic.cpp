#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "m128_utils.h"
#include "generic.h"
#include "m_reg.h"

// t2::t1::t0 -= s2::s1::s0
static uint64_t generic_sub(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
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
static uint64_t generic_add(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
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
static void generic_shr(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s)
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
static void generic_shl(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s)
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
static uint64_t generic_cmp(uint64_t t0, uint64_t t1, uint64_t t2, uint64_t s0, uint64_t s1, uint64_t s2)
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
static uint64_t generic_cmp_sub(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2, borrow;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	borrow = generic_sub(&r0, &r1, &r2, s0, s1, s2);
	if (!borrow) {
		*t0 = r0;
		*t1 = r1;
		*t2 = r2;
	}
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

// binary or
// t2::t1::t0 |= s2::s1:s0 
static void generic_or(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	*t0 = r0 | s0;
	*t1 = r1 | s1;
	*t2 = r2 | s2;
}

static void generic_reduce(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t barrett_lo, uint64_t barrett_hi, uint64_t mod_lo,
			   uint64_t mod_hi)
{
	// process Barrett coefficient
	uint128_t ee_hi, ee_lo, cc, ss0, ss1, ss2;
	uint64_t e_hi, e_lo, s0, s1, s2, borrow;
	ee_hi = (uint128_t) barrett_hi **t2;
	ee_lo = (uint128_t) barrett_hi **t1;
	ee_hi += (ee_lo >> 64);
	ee_lo = (uint128_t) barrett_lo **t2;
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
	borrow = generic_sub(t0, t1, t2, s0, s1, s2);
#if PARANOID
	assert(borrow == 0);
#endif
	while (*t2) {
		generic_sub(t0, t1, t2, mod_lo, mod_hi, 0);
	}
}

static void generic_square(uint128_t v, uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t * t3)
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

// r1::r0 = 2^192 / mod_hi::mod_lo
static void unaligned_reciprocal(uint64_t * r_lo, uint64_t * r_hi, uint64_t mod_lo, uint64_t mod_hi)
{
	uint64_t r0, r1, r2, t0, t1, t2;
	uint64_t k;
	uint128_t q = 0;

	t0 = (uint64_t) - 1;
	t1 = (uint64_t) - 1;
	t2 = (uint64_t) - 1;

	r0 = (uint64_t) 0;
	r1 = (uint64_t) mod_lo;
	r2 = (uint64_t) mod_hi;

	// check division by 0 and even numbers
#if PARANOID
	assert(mod_hi != 0);
	assert((mod_lo & 1) == 1);
#endif

	// align modulus msb on bit 191
	if (r2) {
		k = my_clz64(r2);
	} else if (r1) {
		k = 64 + my_clz64(r1);
	} else {
		k = 128 + my_clz64(r0);
	}
	generic_shl(&r0, &r1, &r2, k);

	// subtract and shift as a double_add implementation
	while (generic_cmp(r0, r1, r2, mod_lo, mod_hi, 0) == 0) {
		q += q;		// quotient double
		if (generic_cmp_sub(&t0, &t1, &t2, r0, r1, r2) == 0) {
			q += 1;	// quotient add
		}
		generic_shr(&r0, &r1, &r2, 1);
	}

	// now
	// r < mod
	// t = (2^192-1) % mod
	// q = 2^192 / mod

	*r_lo = (uint64_t) q;
	*r_hi = (uint64_t) (q >> 64);
}

static uint128_t generic_square_reduce(uint128_t v, uint64_t mod_lo, uint64_t mod_hi)
{
#if PARANOID
	assert(mod_hi > 0);	// at least bit modulus
	assert(mod_lo & 1);	// odd modulus
#endif
	// square first, reduce after
	uint64_t t0, t1, t2, t3, borrow;
	generic_square(v, &t0, &t1, &t2, &t3);

	// Barrett coefficient 2^192 / n
	uint64_t barrett_lo;
	uint64_t barrett_hi;
	unaligned_reciprocal(&barrett_lo, &barrett_hi, mod_lo, mod_hi);

	// modular reduction
	// printf("r2  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 
	generic_reduce(&t1, &t2, &t3, barrett_lo, barrett_hi, mod_lo, mod_hi);
	// printf("r1  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 
	generic_reduce(&t0, &t1, &t2, barrett_lo, barrett_hi, mod_lo, mod_hi);
	// printf("r0  0x%16.16lx %16.16lx %16.16lx %16.16lx\n", t3, t2, t1, t0); 

	// subtract the modulus until the result is less than modulus
	// this loop should run 0, 1 or 2 times, not much more (t < 2 * modulus for sure)
	do {
		// subtract. when a borrow occurs, it was too late.
		borrow = generic_cmp_sub(&t0, &t1, &t2, mod_lo, mod_hi, 0);
	}
	while (!borrow);

	v = ((uint128_t) t1 << 64) + t0;
	return v;
}

static uint128_t genericModSqr(uint128_t v, uint64_t mod_lo, uint64_t mod_hi)
{
	// 1 modulus step
	if ((v >> 64 == 0) && (mod_hi == 0)) {
		uint64_t t_lo = (uint64_t) v;
		v = (uint128_t) t_lo *t_lo;	// < 128 bits
		v %= mod_lo;	// < 64 bits
		return v;
	}
	// 2 modulus steps
	if ((v >> 79) == 0 && (mod_hi < 4)) {
		uint128_t n = ((uint128_t) mod_hi << 64) + mod_lo;
		// - modulus is up to 66 bits
		// - input number is up to 79 bits
		// - result is less than the modulus
		// split, square and reduce without overflow
		//
		// ((((t_hi^2 << 32) + 2*t_lo*t_hi) << 32) + t_lo * t_lo) % n;
		//
		uint64_t t_lo = v & ((1ull << 32) - 1);	// 32 bits
		uint64_t t_hi = v >> 32;	// 47 bits
		v = (uint128_t) t_hi *t_hi;	// 94 bits
		v <<= 32;	// 126 bits
		v += ((uint128_t) t_lo * t_hi) * 2;	// (32+47+1 = 80, 126) --> 127 bits
		v %= n;		// 66 bits
		v <<= 32;	// 98 bits
		v += (uint128_t) t_lo *t_lo;	// (98, 32+32 = 64) --> 99 bits;
		v %= n;		// 66 bits
		return v;
	}
	// 3 modulus steps
	if ((v >> 104) == 0 && (mod_hi >> 23) == 0) {
		// assume n < 87 bits
		// assume v < 104 bits
		// split, square and reduce without overflow
		// return v < n
		//
		uint128_t n = ((uint128_t) mod_hi << 64) + mod_lo;
		uint64_t t_lo = v & ((1ull << 40) - 1);	// 40 bits
		uint64_t t_hi = v >> 40;	// 104 - 40 = 64 bits
		v = (uint128_t) t_hi *t_hi;	// 128 bits
		v %= n;		// 87 bits
		v <<= 40;	// 127 bits;
		v += ((uint128_t) t_lo * t_hi) * 2;	// (127)(105) = 128 bits
		v %= n;		// 87 bits
		v <<= 40;	// 127 bits
		v += (uint128_t) t_lo *t_lo;	// (127)(96) = 128 bits;
		v %= n;		// 87 bits
		return v;
	}
	// 5 modulus steps
	if ((v >> 112) == 0 && (mod_hi >> 39) == 0) {
		// assume n < 103 bits
		// assume v < 112 bits
		// split, square and reduce without overflow
		// return v < n
		//
		// ((((t_hi^2 << 48) + 2*t_lo*t_hi) << 48) + t_lo*t_lo) % n;
		//
		uint128_t n = ((uint128_t) mod_hi << 64) + mod_lo;
		uint64_t t_lo = v & ((1ull << 48) - 1);	// 48 bits
		uint64_t t_hi = v >> 48;	// 112 - 48 = 64 bits
		v = (uint128_t) t_hi *t_hi;	// 128 bits
		v %= n;		// 103 bits
		v <<= 24;	// 127 bits;
		v %= n;		// 103 bits
		v <<= 24;	// 127 bits;
		v += ((uint128_t) t_lo * t_hi) * 2;	// (127)(113) = 128 bits
		v %= n;		// 103 bits
		v <<= 24;	// 127 bits
		v %= n;		// 103 bits
		v <<= 24;	// 127 bits
		v += (uint128_t) t_lo *t_lo;	// (127)(96) = 128 bits;
		v %= n;		// 103 bits
		return v;
	}
	// 9 modulus steps
	if ((v >> 120) == 0 && (mod_hi >> 49) == 0) {
		// assume n < 113 bits
		// assume v < 120 bits
		// split, square and reduce without overflow
		// return v < n
		//
		// ((((t_hi^2 << 48) + 2*t_lo*t_hi) << 48) + t_lo*t_lo) % n;
		//
		uint128_t n = ((uint128_t) mod_hi << 64) + mod_lo;
		uint64_t t_lo = v & ((1ull << 56) - 1);	// 56 bits
		uint64_t t_hi = v >> 56;	// 120 - 56 = 64 bits
		v = (uint128_t) t_hi *t_hi;	// 128 bits
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v += ((uint128_t) t_lo * t_hi) * 2;	// (127)(121) = 128 bits
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v %= n;		// 113 bits
		v <<= 14;	// 127 bits;
		v += (uint128_t) t_lo *t_lo;	// (127)(112) = 128 bits;
		v %= n;		// 113 bits
		return v;
	}

	v = generic_square_reduce(v, mod_lo, mod_hi);
	return v;
}

// Fermat test with Euler criterion.
bool genericFermatTest(uint64_t n_lo, uint64_t n_hi)
{
	uint128_t n = ((uint128_t) n_hi << 64) + n_lo;

#if PARANOID
	// assume n odd
	assert((n & 1) == 1);
#endif

	uint128_t result = 2;
	int bit = 128;
	// find  the msb
	while (bit > 1) {
		bit--;
		if ((n >> bit) & 1) {
			break;
		}
	}
	while (bit > 1) {
		bit--;
		// square
		result = genericModSqr(result, n_lo, n_hi);

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
			if (result >= n) {
				result -= n;
			}
		}
	}
	// Euler criterion
	uint64_t legendre = ((n_lo >> 1) ^ (n_lo >> 2)) & 1;	// shortcut calculation of legendre symbol
	return result % n == (legendre ? n - 1 : 1);
}
