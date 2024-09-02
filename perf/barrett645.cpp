
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>

#include "m_reg.h"

typedef unsigned __int128 uint128_t;

static uint8_t add64(uint64_t * t0, uint64_t * t1, uint64_t s0, uint64_t s1)
{
	uint8_t c;
	uint64_t r0, r1, r2;

	r0 = *t0;
	r1 = *t1;
	c = my_adc64(0, r0, s0, t0);
	c = my_adc64(c, r1, s1, t1);
	return c;		// carry is 1 when t + overflow
}

static uint8_t sub64(uint64_t * t0, uint64_t * t1, uint64_t s0, uint64_t s1)
{
	uint8_t borrow;
	uint64_t r0, r1, r2;

	r0 = *t0;
	r1 = *t1;
	borrow = my_sbb64(0, r0, s0, t0);
	borrow = my_sbb64(borrow, r1, s1, t1);
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

static uint8_t sub64(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint8_t borrow;
	uint64_t r0, r1, r2;

	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	borrow = my_sbb64(0, r0, s0, t0);
	borrow = my_sbb64(borrow, r1, s1, t1);
	borrow = my_sbb64(borrow, r2, s2, t2);
	return borrow;		// borrow is 1 when t < s, borrow is 0 when t >= s
}

static uint8_t cmp_sub64(uint64_t * t0, uint64_t * t1, uint64_t s0, uint64_t s1)
{
	uint64_t r0, r1, r2;
	uint8_t borrow;
	r0 = *t0;
	r1 = *t1;
	borrow = sub64(&r0, &r1, s0, s1);
	*t0 = borrow ? *t0 : r0;
	*t1 = borrow ? *t1 : r1;
	return borrow;
}

static uint8_t cmp_sub64(uint64_t * t0, uint64_t * t1, uint64_t * t2, uint64_t s0, uint64_t s1, uint64_t s2)
{
	uint64_t r0, r1, r2;
	uint8_t borrow;
	r0 = *t0;
	r1 = *t1;
	r2 = *t2;
	borrow = sub64(&r0, &r1, &r2, s0, s1, s2);
	*t0 = borrow ? *t0 : r0;
	*t1 = borrow ? *t1 : r1;
	*t2 = borrow ? *t2 : r2;
	return borrow;
}

// get exact values for barrett = 2^128 / mod and mod128 = 2^128 % mod
static void getBarrett(uint64_t mod_lo, uint64_t * barrett, uint64_t * mod128_lo, uint64_t * mod128_hi)
{
	if (0) {
		// too slow by a few dozen cycles.
		uint128_t mod = ((uint128_t) 1 << 64) + mod_lo;
		// q = 2^128 / mod computed as (2^128-1) / mod  (always the same result)
		uint64_t q;
		uint128_t t = (uint128_t) - 1;
		t /= mod;
		q = (uint64_t) t;
		*barrett = q;
		// 2^128 % mod s computed as 2^128 - mod * 
		t = -(mod * q);
		*mod128_lo = (uint64_t) t;
		*mod128_hi = (uint64_t) (t >> 64);
		return;
	}
	if (mod_lo == 1) {
		// 2^64 + 1 creates an exception in the code below
		*barrett = (uint64_t) - 1;
		*mod128_lo = 1;
		*mod128_hi = 0;
		return;
	}
	// Taken from Hacker's Delight;
	uint128_t c1, c2;
	uint64_t q, mod_hi, r_lo, t0, t1;
	bool c3, c4;

	uint128_t mod = (((uint128_t) 1 << 64) + mod_lo) << 63;
	mod_hi = (uint64_t) (mod >> 64);

	q = my_div64(1ull << 63, 0, mod_hi, &r_lo);
	c1 = (uint128_t) q << 63;
	c2 = (uint128_t) r_lo << 64;
	c3 = c1 > c2;
	c4 = (c1 - c2) > mod;
	q -= c3 ? (c4 ? 2 : 1) : 0;

	t0 = my_mul64(mod_lo, q, &t1);
	t1 += q;
	// 2^128 % (1<<64 + mod_lo)
	uint8_t b = my_sbb64(0, 0, t0, mod128_lo);
	my_sbb64(b, 0, t1, mod128_hi);
	// 2^128 / (1<<64 + mod_lo)
	*barrett = (uint64_t) q;
	return;
}

static void fastSquare(uint64_t * n_lo, uint64_t * n_hi)
{
	// create the correct 128 least significant bits of the squaring
	// result_lo <= 2^63-1 and result_hi is either 0 or 1
	// then this part of the square is less then 2^128
	//
	// the missing part is result_hi*result_hi, which is either 2^128 or 0.
	//
	// This requires the number to square to be less than the modulus
	uint64_t square_lo, square_hi;
	uint64_t r_lo = *n_lo;
	uint64_t r_hi = *n_hi;
	// r_lo * r_lo, r_lo < 2^63
	square_lo = my_mul64(r_lo, r_lo, &square_hi);
	// cross-product 2*r_lo * r_hi, r_hi is either 0 or 1
	square_hi += (r_lo << 1) & (0ull - r_hi);
	// remember the 1 msbit of the squaring over 64.5 bits
	*n_lo = square_lo;
	*n_hi = square_hi;
}

static void fastMod(uint64_t * n_lo, uint64_t * n_hi, uint64_t mod_lo, uint64_t barrett)
{
	// Barrett reduction
	// r = ((n_hi * barrett) >> 64) * mod
	uint64_t e, r0, r1;
	uint64_t n0 = *n_lo;
	uint64_t n1 = *n_hi;
	my_mul64(n1, barrett, &e);
	r0 = my_mul64(e, mod_lo, &r1);
	r1 += e;
	// Barrett reduction
	// n -= r
	sub64(&n0, &n1, r0, r1);
	// conditionally subtract the modulus 2 time
	cmp_sub64(&n0, &n1, mod_lo, 1);
	cmp_sub64(&n0, &n1, mod_lo, 1);
	*n_lo = n0;
	*n_hi = n1;
}

bool fermatTest645(uint64_t n_lo)
{
	// n must be < 2^64.5 ~ 26087635650665564424 = 2.6087e19
	uint64_t bit, result_lo, result_hi, overflow645;
	uint64_t barrett, mod128_lo, mod128_hi;

	getBarrett(n_lo, &barrett, &mod128_lo, &mod128_hi);

	bit = 63;
	uint64_t s, msbs = 2;

	// At least 5 iterations without multiplication, 
	// and without modular reduction
	while (bit) {
		bit -= 1;
		s = msbs;	// save the current bit position
		msbs *= 2;
		msbs += (n_lo >> bit) & 1;
		if (msbs >= 63)
			break;
	}
	// rollback the last iteration
	bit += 1;
	// s will be strictly less than 63
	// do the exponentiation 2^s mod n
	result_lo = 1ull << s;
	result_hi = 0;
	// result_lo is now a 62 bit number, i.e. result < n (the modulus)

	while (bit > 1) {
		bit -= 1;

		// start squaring 
		overflow645 = result_hi;
		fastSquare(&result_lo, &result_hi);
		// reduce
		fastMod(&result_lo, &result_hi, n_lo, barrett);

		if (overflow645) {
			// finish squaring
			add64(&result_lo, &result_hi, mod128_lo, mod128_hi);
			// reduce
			cmp_sub64(&result_lo, &result_hi, n_lo, 1);
		}
		if ((n_lo >> bit) & 1) {
			// double
			my_shld64(&result_hi, &result_lo, 1);
			// reduce
			cmp_sub64(&result_lo, &result_hi, n_lo, 1);
		}
	}

	if (((n_lo >> 1) ^ (n_lo >> 2)) & 1) {
		return result_lo == n_lo - 1 && result_hi == 1;
	} else {
		return result_lo == 1 && result_hi == 0;
	}
}

void self_tests(void)
{
	uint128_t m, v, t;
	uint64_t barrett, mod128_lo, mod128_hi, v_lo, v_hi, overflow, t_lo, t_hi;

	// min modulus for barrett
	m = (uint128_t) 1 << 64;
	m += 1;
	getBarrett((uint64_t) m, &barrett, &mod128_lo, &mod128_hi);
	// printf("%16.16lx %16.16lx %16.16lx\n", barrett, mod128_lo, mod128_hi);
	//
	assert(barrett == (uint64_t) - 1);
	assert(mod128_lo == 1);
	assert(mod128_hi == 0);

	// max modulus for barrett
	m = (uint128_t) 1 << 64;
	m += (1ull << 63) - 1;;
	getBarrett((int64_t) m, &barrett, &mod128_lo, &mod128_hi);
	// printf("%16.16lx %16.16lx %16.16lx\n", barrett, mod128_lo, mod128_hi);
	assert(barrett == 0xaaaaaaaaaaaaaaabull);
	assert(mod128_lo == 0x2aaaaaaaaaaaaaabull);
	assert(mod128_hi == 0);

	// min modulus for square
	m = (uint128_t) 1 << 64;
	m += 1;
	getBarrett((uint64_t) m, &barrett, &mod128_lo, &mod128_hi);
	// max_value to square
	v_hi = 1;
	v_lo = 0x6a09e667f3bcc907ull;
	v = ((uint128_t) v_hi << 64) + v_lo;
	overflow = v_hi;
	fastSquare(&v_lo, &v_hi);
	// printf("%16.16lx %16.16lx\n", v_lo, v_hi);
	assert(v_lo == 0x31b0487d2a23fe31ull);
	assert(v_hi == 0xfffffffffffffffbull);
	fastMod(&v_lo, &v_hi, m, barrett);
	// printf("%16.16lx %16.16lx\n", v_lo, v_hi);
	assert(v_lo == 0x31b0487d2a23fe37ull);
	assert(v_hi == 0);
	if (overflow) {
		add64(&v_lo, &v_hi, mod128_lo, mod128_hi);
		// reduce
		cmp_sub64(&v_lo, &v_hi, m, 1);
	}
	// printf("%16.16lx %16.16lx\n", v_lo, v_hi);
	assert(v_lo == 0x31b0487d2a23fe38ull);
	assert(v_hi == 0);

	// v^2 mod m, slow version
	t = v >> 32;
	v = v & 0xffffffff;
	t = ((((((t * t) << 32) + (t * v * 2)) % m) << 32) + v * v) % m;
	t_lo = (uint64_t) t;
	t_hi = (uint64_t) (t >> 64);
	// printf("%16.16lx %16.16lx\n", t_lo, t_hi);
	assert(t_lo == v_lo);
	assert(t_hi == v_hi);
}

static uint64_t seed = 0x87654321abcdef01ull;
uint64_t my_random(void)
{
	seed = seed * 137 + 13;
	return (seed >> 13) ^ (seed << 13) ^ seed;
}

int main(int argc, char **argv)
{
	// destroys the optimizer which could unroll and strip of the code
	printf("%s: Fermat-Euler test base 2 with simplified squaring and random modulus.\n", argv[0]);
	if (argc > 1)
		seed = time(NULL);
	self_tests();
	fflush(stdout);
	fflush(stderr);

	volatile uint64_t t = 0;
	uint64_t count = 0;
	uint64_t pseudoprimes = 0;
	uint64_t mod_min = 0x0000000000000001ull;	// included
	uint64_t mod_max = 0x6a09e667f3bcc908ull;	// excluded
	uint64_t mod_range = mod_max - mod_min;
	bool b645 = false;
	for (count = 0; count < 2000000; count++) {
		// make a odd number with 2 top msbits = 0b10
		uint64_t n_lo = ((my_random() % mod_range) + mod_min) | 1;
		t -= my_rdtsc();
		b645 = fermatTest645(n_lo);
		t += my_rdtsc();
		pseudoprimes += b645 ? 1 : 0;
	}
	t /= count;
	printf("%lu pseudoprimes / %lu integers, %lu cycles/fermat.\n", pseudoprimes, count, t);
	return 0;
}
