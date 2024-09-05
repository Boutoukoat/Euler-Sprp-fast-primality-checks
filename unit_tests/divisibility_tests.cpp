#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "divisibility.h"

// random proven primes
static uint64_t mersennes[] = { 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 43112609, 57885161, 74207281, 77232917, 82589933 };

static uint64_t factors[] = { 13, 31, 41, 61, 151, 17, 29, 43, 113, 127, 19, 37, 73, 109, 23, 89 };
static uint64_t n_factors[] = { 49, 53, 59, 67, 101, 103 };

int main(int argc, char **argv)
{
	bool verbose = false;
	uint64_t check_count = 0;

	for (int i = 1; i < argc; i++)
	{
		if (!strcmp(argv[i], "-v"))
		{
			verbose = true;
		}
	}

	// The divisibility code is full of constants from nowhere (described in compiler tech papers),
	// better verify there is no corruption of one of the constants
	
	for (unsigned i = 0; i < sizeof(mersennes)/sizeof(mersennes[0]); i++)
	{
		// get 128 bit number with large factors
		uint64_t m = mersennes[i];
		m *= m;
		uint128_t n = m;
		n *= n;

		for (unsigned j = 0; j < sizeof(factors)/sizeof(factors[0]); j++)
		{
			uint64_t t_lo, t_hi;
			bool bb;
			uint128_t t = n * factors[j];
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == false);
			check_count += 1;
			t = n;
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == true);
			check_count += 1;
			t = n * 3 * 5 * 7 * 11;
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == true);
			check_count += 1;
			t = n * factors[j] * factors[j];
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == false);
			check_count += 1;
			t = (uint128_t)1 << 127;
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == true);
			check_count += 1;
		}
		for (unsigned j = 0; j < sizeof(n_factors)/sizeof(n_factors[0]); j++)
		{
			uint64_t t_lo, t_hi;
			bool bb;
			uint128_t t = n * n_factors[j];
			t_lo = (uint64_t)t;
			t_hi = (uint64_t)(t >> 64);
			bb = divisibility_sieve(t_lo, t_hi);
			assert(bb == true);
			check_count += 1;
		}
		if (verbose)
		{
			printf("%lu tests passed\n", check_count);
		}
	}

	printf("%lu tests passed\n", check_count);
	printf("All tests passed\n");
	return 0;
}


