#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "m_reg.h"
#include "m128_utils.h"
#include "slow.h"
#include "generic.h"
#include "barrett.h"
#include "montgomery.h"

static uint128_t my_random(void)
{
	// based on congruential generator, period = 2^128
	static uint128_t seed = ((uint128_t) 0x123456789ull << 92) + ((uint128_t) 0xabcdef << 36) + 0x967654321ull;
	seed = seed * 137 + 13;
	// shuffle
	uint128_t x = seed ^ (seed >> 17) ^ (seed << 17);
	return x;
}

int main(int argc, char **argv)
{
	// variables declared volatile to prevent the optimizer
	// to play with them and foul the timing
	volatile uint64_t tgeneric = 0;
	volatile uint64_t tslow = 0;
	volatile uint64_t tbarrett = 0;
	volatile bool b = 0;
	volatile uint64_t t = my_rdtsc();
	bool bbarrett = false, bgeneric = true, bslow = true;
	bool csv = false;
	uint64_t start_bit = 3;
	uint64_t end_bit = 128;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-start")) {
			start_bit = strtol(argv[++i], 0, 0);
			continue;
		}
		if (!strcmp(argv[i], "-end")) {
			end_bit = strtol(argv[++i], 0, 0);
			continue;
		}
		if (!strcmp(argv[i], "-csv")) {
			csv = true;
			continue;
		}
		if (!strcmp(argv[i], "-nocsv")) {
			csv = false;
			continue;
		}
		printf("%s -start xxx - end xxxx [-csv | -nocsv]\n", argv[0]);
		exit(1);
	}

	if (!csv) {
		printf("Testing bit range from %ld (included) to %ld (included)\n", start_bit, end_bit);
	}

	if (start_bit < 3) {
		printf("prime numbers from 0 to 3 are very hard to detect\n");
		exit(1);
	}
	if (end_bit > 128) {
		printf("prime numbers larger than 128 bits are very hard to detect\n");
		exit(1);
	}
	if (end_bit < start_bit) {
		printf("start of the range > end of the range\n");
		exit(1);
	}

	fflush(stdout);
	fflush(stderr);

	if (!csv) {
#if defined(__x86_64__)
		printf("Compiled for X86_64 target\n");
#elif defined (__aarch64__)
		printf("Compiled for ARM64 target\n");
#else
		printf("Compiled for an unknown target\n");
#endif

#ifdef INTEL_COMPILER
		printf("Compiled with ICC\n");
#elif defined(__GNUC__)
#if defined(__clang__)
		printf("Compiled with CLANG\n");
#else
		printf("Compiled with GCC\n");
#endif
#else
		printf("Unknown compiler\n");
#endif

#if defined(__OPTIMIZATION_LEVEL__)
		printf("Compiled with optimization level -O%d\n", __OPTIMIZATION_LEVEL__);
#endif

#if PARANOID
		printf("PARANOID excessive checks enabled (slow)\n");
#endif
	} else {
		printf("\"bits\";\"slow\";\"generic\";\"barrett\";\n");
	}

	for (uint64_t bits = start_bit; bits <= end_bit; bits++) {
		if (!csv) {
			printf("----- %3lu bits:\n", bits);
		}
		fflush(stdout);
		fflush(stderr);

		uint128_t x;
		do {
			x = my_random();
			uint128_t maskAnd = ((uint128_t) 1 << (bits - 1)) - 1;	// clear msbits
			uint128_t maskOr = ((uint128_t) 1 << (bits - 1)) | ((uint128_t) 1 << (bits / 2));	// force msb, force another bit
			x &= maskAnd;
			x |= maskOr;
			x /= 6;
			x *= 6;	// now a multiple of 6
			x += 1;	// number like 6*k + 1
		}
		while (x >> (bits - 1) != 1);

		// printf("start x "); my_printf(x); printf("\n");
		uint128_t xend = (bits == 128) ? (uint128_t) - 1 : ((uint128_t) 1 << bits);
		uint128_t xstart = x;
		uint64_t x_lo, x_hi;
		uint64_t xinc = 4;
		uint64_t count = 0;
		uint64_t countPrime = 0;

		tbarrett = 0;
		tgeneric = 0;
		tslow = 0;
		count = 0;

		// Instruction Cache warmup 
		x_lo = (uint64_t) x;
		x_hi = (uint64_t) (x >> 64);
		t -= my_rdtsc();
		b ^= genericFermatTest(x_lo, 0x0);
		b |= barrettFermatTest(x_lo, 0x0);
		b &= my_slowEuler2(x);
		b ^= genericFermatTest(x_lo, x_hi);
		b |= barrettFermatTest(x_lo, x_hi);
		b &= my_slowEuler2(x);
		t += my_rdtsc();
		b &= (t == 1);	// some stupid dependency to prevent the compiler to strip out the code

		// minimum test dalay
		const uint64_t delay = 3e8;

		// run for at least 1000000000 cycles of testing
		while (tbarrett < delay || tgeneric < delay || tslow < delay) {
			// make sure to test only 'bits' number
			while (x < xend && (tbarrett < delay || tgeneric < delay || tslow < delay)) {
				// test numbers 6^k -1
				x_lo = (uint64_t) x;
				x_hi = (uint64_t) (x >> 64);
				t = my_rdtsc();
				tslow -= t;
				bslow = my_slowEuler2(x);
				b ^= bslow;
				t = my_rdtsc();
				tslow += t;
				tgeneric -= t;
				bgeneric = genericFermatTest(x_lo, x_hi);
				b ^= bgeneric;
				t = my_rdtsc();
				tgeneric += t;
				tbarrett -= t;
				bbarrett = barrettFermatTest(x_lo, x_hi);
				b ^= bbarrett;
				t = my_rdtsc();
				tbarrett += t;

				if (bgeneric != bslow || bbarrett != bslow) {
					// some error occured
					// aoutch.  

					// get a third party consensus
					bool f2slow = my_slowFermat2(x);
					bool f3slow = my_slowFermat3(x);

					// display booleans
					printf("booleans ..... : %d %d %d %d %d\n",
					       (int)bgeneric, (int)bbarrett, (int)bslow, (int)f2slow, (int)f3slow);

					// display the result of many different tests
					printf("modulus ...... : ");
					my_printf(x);
					printf("\n");
					printf("\n");
					printf("Fermat-Euler tests\n");
					printf("barrett ...... : %s\n", bbarrett ? "might be prime" : "composite for sure");
					printf("generic ...... : %s\n", bgeneric ? "might be prime" : "composite for sure");
					printf("slow ......... : %s\n", bslow ? "might be prime" : "composite for sure");
					printf("\n");
					printf("Fermat tests\n");
					printf("base 2 ....... : %s\n", f2slow ? "might be prime" : "composite for sure");
					printf("base 3 ....... : %s\n", f3slow ? "might be prime" : "composite for sure");
					printf("\n");
					printf("Fail\n");
					assert(bslow == bgeneric);
					abort();
				}
				x += xinc;
				xinc = 6 - xinc;
				count += 1;
				countPrime += bgeneric;
			}
			// x might not be a 'bits' number, restart if needed from the start of the range
			x = xstart;
			xinc = 4;
		}

		// display cycles and averages
		double avgslow = (double)tslow / count;
		double avggeneric = (double)tgeneric / count;
		double avgbarrett = (double)tbarrett / count;
		if (csv) {
			printf("%ld;%.1f;%.1f;%.1f;\n", bits, avgslow, avggeneric, avgbarrett);
		} else {
			printf("slow .......... : %12.1f cycles / fermat\n", avgslow);
			printf("generic ....... : %12.1f cycles / fermat\n", avggeneric);
			printf("barrett ....... : %12.1f cycles / fermat\n", avgbarrett);
			printf("%ld pseudoprimes, %ld composites\n", countPrime, count - countPrime);
		}
	}

	fflush(stdout);
	fflush(stderr);
	return b;		// add some dependency to prevent 
	// the compiler/optimizer to strip the whole code out
}
