#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "slow.h"

// shift and add
// return value < mod
uint128_t my_slowReduce(uint128_t x, uint128_t mod)
{
	// reduce for all possible inputs > 0  and < 128 bits
	if (x < mod) {
		return x;
	}
	x -= mod;
	if (x < mod) {
		return x;
	}
	// remove some multiples of the modulus
	if (mod >> 124 == 0 && x >> 127 == 0) {
		uint128_t u = mod;
		while (u < x) {
			u <<= 1;
		}
		while (u >= mod) {
			if (x >= u) {
				x -= u;
			}
			u >>= 1;
		}
		return x;
	}
	// remove iteratively the modulus
	do {
		x -= mod;
	}
	while (x <= mod);
	return x;
}

// modular addition for mod < 128 bits, a < mod and b < mod
// return value < mod
uint128_t my_slowModAdd(uint128_t a, uint128_t b, uint128_t mod)
{
	uint128_t t = a + b;
	if (t < a) {
		// wrap-around occured , this could occur if mod is 127 bits.
		// input condition : a <= n-1, b <= n-1  then  (a+b) <= 2*n-2
		uint128_t u = t;
		t = t - mod;	// subtract the modulus, ignore carry and borrow.
#if PARANOID
		// already checked carry was 1 and now check borrow is 1
		// if this fails, this means a < mod and b < mod are not true and something is wrong
		assert(t >= u);	// another wrap-around must occur
		assert(t < mod);	// reduction is complete because (a+b) - n <= n-2

#endif
		return t;
	}
	// most frequent case
	if (t < mod) {
		return t;
	}
	t -= mod;
	if (t < mod) {
		return t;
	}
	// unreachable since a < mod and b < mod
	return my_slowReduce(t, mod);
}

// shift and add squaring for all possible x, mod inputs, return value < mod
uint128_t my_slowModSqr(uint128_t x, uint128_t mod)
{
	// square and reduce for all possible inputs < 128 bits
	uint128_t result = 0;
	int bit = 128;
	x = my_slowReduce(x, mod);

	// skip quickly the zeroes
	while (bit > 1) {
		bit--;
		if ((x >> bit) & 1) {
			result = x;
			break;
		}
	}

	while (bit--) {
		// double and reduce
		result = my_slowModAdd(result, result, mod);

		if ((x >> bit) & 1) {
			// add and reduce
			result = my_slowModAdd(result, x, mod);
		}
	}
	return result;
}

// Fermat test with Euler criterion
bool my_slowEuler2(uint128_t mod)
{
#if PARANOID
	assert((mod & 1) == 1);
#endif
	uint128_t result = 1;
	int bit = 128;
	while (bit > 1) {
		bit--;
		if ((mod >> bit) & 1) {
			// multiply by 2
			result = 2;
			break;
		}
	}

	while (bit > 1) {
		bit--;
		// square
		result = my_slowModSqr(result, mod);
		if ((mod >> bit) & 1) {
			// multiply by 2
			result = my_slowModAdd(result, result, mod);
		}
	}

	// Euler criterion
	uint64_t legendre = ((mod >> 1) ^ (mod >> 2)) & 1;	// shortcut calculation of legendre symbol
	return (result == (legendre ? mod - 1 : 1));
}

// Fermat test base 2 
bool my_slowFermat2(uint128_t mod)
{
#if PARANOID
	assert((mod & 1) == 1);
#endif
	uint128_t result = 1;
	int bit = 128;
	while (bit > 1) {
		bit--;
		if ((mod >> bit) & 1) {
			result = 2;
			break;
		}
	}

	while (bit > 1) {
		bit--;
		// square
		result = my_slowModSqr(result, mod);
		if ((mod >> bit) & 1) {
			// multiply by 2
			result = my_slowModAdd(result, result, mod);
		}
	}
	// square
	result = my_slowModSqr(result, mod);
	return (result == 1);
}

// Fermat test base 3 
bool my_slowFermat3(uint128_t mod)
{
#if PARANOID
	assert((mod & 1) == 1);
#endif
	uint128_t result = 1;
	int bit = 128;
	while (bit > 1) {
		bit--;
		if ((mod >> bit) & 1) {
			result = 3;
			break;
		}
	}

	while (bit > 1) {
		bit--;
		// square
		result = my_slowModSqr(result, mod);
		if ((mod >> bit) & 1) {
			// multiply by 3
			uint128_t t = my_slowModAdd(result, result, mod);
			result = my_slowModAdd(result, t, mod);
		}
	}
	// square
	result = my_slowModSqr(result, mod);
	return (result == 1);
}

// SPRP base 2 test 
bool my_slowSprp2(uint128_t mod)
{
#if PARANOID
	assert((mod & 1) == 1);
#endif
	uint128_t result = 1;
	uint128_t k = 0;
	uint128_t scan = mod - 1;
	while ((scan & 1) == 0) {
		scan >>= 1;
		k += 1;
	}
	int bit = 128;
	while (bit > 1) {
		bit--;
		if ((scan >> bit) & 1) {
			result = 2;
			break;
		}
	}

	while (bit > 0) {
		bit--;
		// square
		result = my_slowModSqr(result, mod);
		if ((scan >> bit) & 1) {
			// multiply by 2
			result = my_slowModAdd(result, result, mod);
		}
	}

	if (result == 1) {
		return true;
	}
	while (k > 1) {
		k -= 1;
		if (result == mod - 1) {
			return true;
		}
		// square
		result = my_slowModSqr(result, mod);

	}
	if (result == mod - 1) {
		return true;
	}
	return false;
}

// SPRP base 3 test 
bool my_slowSprp3(uint128_t mod)
{
#if PARANOID
	assert((mod & 1) == 1);
#endif
	uint128_t result = 1;
	uint128_t k = 0;
	uint128_t scan = mod - 1;
	while ((scan & 1) == 0) {
		scan >>= 1;
		k += 1;
	}
	int bit = 128;
	while (bit > 1) {
		bit--;
		if ((scan >> bit) & 1) {
			result = 3;
			break;
		}
	}

	while (bit > 0) {
		bit--;
		// square
		result = my_slowModSqr(result, mod);
		if ((scan >> bit) & 1) {
			// multiply by 3
			uint128_t t = my_slowModAdd(result, result, mod);
			result = my_slowModAdd(t, result, mod);
		}
	}

	if (result == 1) {
		return true;
	}
	while (k > 1) {
		k -= 1;
		if (result == mod - 1) {
			return true;
		}
		// square
		result = my_slowModSqr(result, mod);

	}
	if (result == mod - 1) {
		return true;
	}
	return false;
}
