#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "m_reg.h"
#include "m64_utils.h"
#include "m128_utils.h"
#include "slow.h"
#include "optimized.h"
#include "avx2_sprp.h"

static uint128_t my_random(void)
{
	// based on linear congruential generator, period = 2^128
	static uint128_t seed = ((uint128_t) 0x123456789ull << 92) + ((uint128_t) 0xabcdef << 36) + 0x987654321ull;
	seed = seed * 137 + 13;
	// shuffle
	uint128_t x = seed ^ (seed >> 17) ^ (seed << 13);
	return x;
}

static bool isprime(uint128_t x)
{
	// good enough for finding primes for performance testing 
	// A rare pseudoprime (without a small factor) would not affect the average timing
	if (x > 5 && x % 5 == 0) return false;
	if (x > 7 && x % 7 == 0) return false;
	if (x > 11 && x % 11 == 0) return false;
	if (x > 13 && x % 13 == 0) return false;
	if (x > 17 && x % 17 == 0) return false;
	if (x > 19 && x % 19 == 0) return false;
	if (x > 23 && x % 23 == 0) return false;
	return my_slowSprp2(x) && my_slowSprp3(x) && my_slowSprp5(x); 
}

static bool avx256_test(uint128_t v)
{
        uint32_t mm = montgomeryInverse32((uint32_t) v);

	if (v >> 32 == 0)
	{
		uint32_t montg_bases[4];
		uint32_t one = montgomery_bases1(montg_bases, v, 4);
		return avx2_sprp1(v, mm, one, montg_bases);
	}

	if (v >> 64 == 0)
	{
		uint64_t montg_bases[4];
		uint64_t one = montgomery_bases2(montg_bases, v, 4);
		return avx2_sprp2(v, mm, one, montg_bases);
	}

	if (v >> 96 == 0)
	{
		uint128_t montg_bases[4];
		uint128_t one = montgomery_bases3(montg_bases, v, 4);
		return avx2_sprp3(v, mm, one, montg_bases);
	}

	else
	{
		uint128_t montg_bases[4];
		uint128_t one = montgomery_bases4(montg_bases, v, 4);
		return avx2_sprp4(v, mm, one, montg_bases);
	}
}

static bool avx512_test(uint128_t v)
{
        uint32_t mm = montgomeryInverse32((uint32_t) v);

	if (v >> 32 == 0)
	{
		uint32_t montg_bases[8];
		uint32_t one = montgomery_bases1(montg_bases, v, 8);
#ifdef __AVX512F__
		return avx512_sprp1(v, mm, one, montg_bases);
#else
		return avx22_sprp1(v, mm, one, montg_bases);
#endif
	}

	if (v >> 64 == 0)
	{
		uint64_t montg_bases[8];
		uint64_t one = montgomery_bases2(montg_bases, v, 8);
#ifdef __AVX512F__
		return avx512_sprp2(v, mm, one, montg_bases);
#else
		return avx22_sprp2(v, mm, one, montg_bases);
#endif
	}

	if (v >> 96 == 0)
	{
		uint128_t montg_bases[8];
		uint128_t one = montgomery_bases3(montg_bases, v, 8);
#ifdef __AVX512F__
		return avx512_sprp3(v, mm, one, montg_bases);
#else
		return avx22_sprp3(v, mm, one, montg_bases);
#endif
	}
	else
	{
		uint128_t montg_bases[8];
		uint128_t one = montgomery_bases4(montg_bases, v, 8);
#ifdef __AVX512F__
		return avx512_sprp4(v, mm, one, montg_bases);
#else
		return avx22_sprp4(v, mm, one, montg_bases);
#endif
	}
}

int main(int argc, char **argv)
{
	// variables declared volatile to prevent the avx512imizer
	// to play with them and foul the timing
	volatile uint64_t tavx256 = 0;
	volatile uint64_t topt = 0;
	volatile uint64_t tavx512 = 0;
	volatile bool b = 0;
	volatile uint64_t t = my_rdtsc();
	bool bavx256 = false, bopt = false, bavx512 = false;
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
		printf("Compiled with avx512imization level -O%d\n", __OPTIMIZATION_LEVEL__);
#endif

#if PARANOID
		printf("PARANOID excessive checks enabled (opt)\n");
#endif
#if defined(__AVX2__)
		printf("Compiled for AVX-256 target\n");
#endif
#if defined(__AVX512F__)
		printf("Compiled for AVX-512 target\n");
#endif
	} else {
#if defined(__AVX512F__)
		printf("\"bits\";\"uint64_t\";\"avx2_epu32\";\"avx512_epu32\";\"avx512_ifma\";\n");
#else
		printf("\"bits\";\"uint64_t\";\"avx2_epu32\";\"avx2_stitched\";\"avx2_ifma\";\n");
#endif
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
		uint64_t countopt = 0;
		uint64_t count256 = 0;
		uint64_t count512 = 0;
		uint64_t countPrime = 0;

		topt = 0;
		tavx256 = 0;
		tavx512 = 0;

		// Instruction Cache warmup 
		x_lo = (uint64_t) x;
		x_hi = (uint64_t) (x >> 64);
		t -= my_rdtsc();
		b ^= optimizedSprpTest(x_lo, x_hi);
		b |= avx2SprpTest(x_lo, x_hi);
		b &= avxSprpTest(x_lo, x_hi);
		b ^= optimizedSprpTest(x_lo, x_hi);
		b |= avx2SprpTest(x_lo, x_hi);
		b &= avxSprpTest(x_lo, x_hi);
		b ^= optimizedSprpTest(x_lo, x_hi);
		b |= avx2SprpTest(x_lo, x_hi);
		b &= avxSprpTest(x_lo, x_hi);
		t += my_rdtsc();
		b &= (t == 1);	// some stupid dependency to prevent the compiler to strip out the code

		// minimum test dalay
		const uint64_t delay = 1e7;

		// run for at least 1000000000 cycles of testing
		while (topt < delay || tavx512 < delay || tavx256 < delay) {
			// make sure to test only 'bits' number
			while (x < xend && (topt < delay || tavx512 < delay || tavx256 < delay)) {
				// test numbers 6^k -1
				bool prime_only = isprime(x);
				if (prime_only)
				{
				x_lo = (uint64_t) x;
				x_hi = (uint64_t) (x >> 64);
				t = my_rdtsc();
				topt -= t;
				asm volatile("#");
				bopt = optimizedSprpTest(x_lo, x_hi);  // uint128_t
				asm volatile("#");
				b ^= bopt;
				t = my_rdtsc();
				topt += t;
				tavx512 -= t;
				asm volatile("#");
				bavx512 = avx512_test(x);    // avx512 or 2 x avx256
				asm volatile("#");
				b |= bavx512;
				t = my_rdtsc();
				tavx512 += t;
				tavx256 -= t;
				asm volatile("#");
				bavx256 = avx256_test(x);     // avx256 only
				asm volatile("#");
				b -= bavx256;
				t = my_rdtsc();
				tavx256 += t;

				if (bavx256 != bavx512) {
					// some error occured
					// aoutch.  

					// get a third party consensus
					bool f2opt = my_slowFermat2(x);
					bool f3opt = my_slowFermat3(x);
					bool bopt2 = my_slowSprp2(x);
					bool bopt3 = my_slowSprp3(x);

					// display booleans for debug
					printf("booleans ..... : %d %d %d %d %d %d\n", (int)bopt3, (int)bopt2, (int)bavx256,
					       (int)bavx512, (int)f2opt, (int)f3opt);

					// display the result of many different tests
					printf("modulus ...... : ");
					my_printf(x);
					printf("\n");
					printf("\n");
					printf("Sprp tests\n");
					printf("avx512 ....... : %s\n", bavx512 ? "might be prime" : "composite for sure");
					printf("avx256 ....... : %s\n", bavx256 ? "might be prime" : "composite for sure");
					printf("\n");
					printf("Slow tests\n");
					printf("base 2........ : %s\n", bopt2 ? "might be prime" : "composite for sure");
					printf("base 3 ....... : %s\n", bopt3 ? "might be prime" : "composite for sure");
					printf("\n");
					printf("Fermat tests\n");
					printf("base 2 ....... : %s\n", f2opt ? "might be prime" : "composite for sure");
					printf("base 3 ....... : %s\n", f3opt ? "might be prime" : "composite for sure");
					printf("\n");
					printf("Fail\n");
					assert(bavx512 == bavx256);
					abort();
				}
				countopt += 1;  // uint128_t
				count256 += 4;  // avx256 
				count512 += 8;  // avx512
				countPrime += bavx512;
				}
				x += xinc;
				xinc = 6 - xinc;
			}
			// x might not be a 'bits' number, restart if needed from the start of the range
			x = xstart;
			xinc = 4;
		}

		// display cycles and averages
		double avgopt = (double)topt / countopt;
		double avgavx256 = (double)tavx256 / count256;
		double avgavx512 = (double)tavx512 / count512;
		if (csv) {
			printf("%ld;%.1f;%.1f;%.1f;\n", bits, avgopt, avgavx256, avgavx512);
		} else {
			printf("uint64_t ........... : %12.1f cycles / sprp\n", avgopt);
			printf("avx256 ............. : %12.1f cycles / sprp\n", avgavx256);
			printf("avx512 ............. : %12.1f cycles / sprp\n", avgavx512);
			printf("%ld pseudoprimes, %ld composites\n", countPrime, countopt - countPrime);
		}
	}

	fflush(stdout);
	fflush(stderr);
	return b;		// add some dependency to prevent 
	// the compiler/avx512imizer to strip the whole code out
}
