#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <x86intrin.h>

#if PARANOID
#include "montgomery.h"
#endif

#include "optimized.h"

#include "avx2_sprp.h"

// utilities to convert strings from/to 128 bits numbers
#include "m128_utils.h"

// divisibility sieve internals
#include "m_reg.h"

// use a simplistic wheel sieve mod 2,3,5,7,11
#include "gap_check.wheel"

// ------------------------------------------------------------------------
// fast divisibility by constants 13, 31, 41, 61 , the barrett way
// ------------------------------------------------------------------------
static bool divisibility_sieve(uint64_t t_lo, uint64_t t_hi)
{
	uint64_t v_lo, v_hi;
        // Use the factors of mersenne number 2^60-1 because the modular reduction is simple
        // The factors of 2^60 - 1 are 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321

        // first, reduce mod (2^60-1)
	v_lo = t_lo;
	v_hi = t_hi;
        uint64_t mask60 = (1ull << 60) - 1;
        uint64_t mod60 = v_lo & mask60;
        my_shrd64(&v_hi, &v_lo, 60);
        mod60 += v_lo & mask60;
        my_shrd64(&v_hi, &v_lo, 60);
        mod60 += v_lo;
        // mod60 is now a 62 bit number with the same factors than the input number v.
        // i.e. gcd(2^60-1, mod60) == gcd(2^60-1, v)

        // since factors 3, 2, 5, 11 are skipped, let's start with 13
	if (mod60 % 13 == 0) return false;  // composite for sure
        if (mod60 % 31 == 0) return false;  // composite for sure
        if (mod60 % 41 == 0) return false;  // composite for sure
        if (mod60 % 61 == 0) return false;  // composite for sure

        // first, reduce mod (2^56-1)
	// 2^56 - 1 = 3 * 5 * 17 * 29 * 43 * 113 * 127 * 15790321
	v_lo = t_lo;
	v_hi = t_hi;
        uint64_t mask56 = (1ull << 56) - 1;
        uint64_t mod56 = v_lo & mask56;
        my_shrd64(&v_hi, &v_lo, 56);
        mod56 += v_lo & mask56;
        my_shrd64(&v_hi, &v_lo, 56);
        mod56 += v_lo;
        // mod56 is now a 58 bit number with the same factors than the input number v.
	if (mod56 % 17 == 0) return false;  // composite for sure
        if (mod56 % 29 == 0) return false;  // composite for sure
        if (mod56 % 43 == 0) return false;  // composite for sure

        return true;            // might be prime
}

// search a prime backwards in range ]lo, hi], from hi to lo
static uint128_t search_backwards(uint128_t lo, uint128_t hi)
{
	// use a simplistic wheel sieve mod 2,3,5,7,11
	// search the first prime candidate using the sieve
	uint32_t t = hi % (sizeof(sieve_backwards) / sizeof(sieve_backwards[0]));
	uint32_t v = sieve_backwards[t - 1];
	if (!sieve_backwards[t]) {
		// composite number, move to the next sieve candidate
		hi -= (t > v) ? t - v : (t + sizeof(sieve_backwards) / sizeof(sieve_backwards[0])) - v;
		t = v;
		v = sieve_backwards[t - 1];
	}
	// verify the prime candidate and move to the next sieve candidate
	while (hi > lo) {
		// verify the prime candidate 
		uint64_t t_lo = (uint64_t) hi;
		uint64_t t_hi = (uint64_t) (hi >> 64);
		bool bd = divisibility_sieve(t_lo, t_hi);
		if (bd) {
			// no trivial factor
			bool bs = optimizedSprpTest(t_lo, t_hi);
			if (bs) {
				// pseudoprime found
				bool ba = avx2SprpTest(t_lo, t_hi);
				if (ba) {
					// prime found
					break;
				}
			}
		}
		// skip composites, move backwards to the next sieve candidate
		hi -= (t > v) ? t - v : (t + sizeof(sieve_backwards) / sizeof(sieve_backwards[0])) - v;
		t = v;
		v = sieve_backwards[t - 1];
	}
	return hi;
}

// search a prime forwards in range [lo, hi[, from lo to hi
static uint128_t search_forwards(uint128_t lo, uint128_t hi)
{
	// use a simplistic wheel sieve mod 2,3,5,7,11
	// search the first prime candidate using the sieve
	uint32_t t = lo % (sizeof(sieve_forwards) / sizeof(sieve_forwards[0]));
	uint32_t v = sieve_forwards[t - 1];
	if (!sieve_forwards[t]) {
		// composite number, move to the next sieve candidate
		lo += (v > t) ? v - t : (v + sizeof(sieve_forwards) / sizeof(sieve_forwards[0])) - t;
		t = v;
		v = sieve_forwards[t - 1];
	}
	// verify the prime candidate and move to the next sieve candidate
	while (lo < hi) {
		// verify the prime candidate 
		uint64_t t_lo = (uint64_t) lo;
		uint64_t t_hi = (uint64_t) (lo >> 64);
		bool bd = divisibility_sieve(t_lo, t_hi);
		if (bd) {
			bool bs = optimizedSprpTest(t_lo, t_hi);
			if (bs) {
				// pseudoprime found
				bool ba = avx2SprpTest(t_lo, t_hi);
				if (ba) {
					// prime found
					break;
				}
			}
		}
		// skip composites, move to the next sieve candidate
		lo += (v > t) ? v - t : (v + sizeof(sieve_forwards) / sizeof(sieve_forwards[0])) - t;
		t = v;
		v = sieve_forwards[t - 1];
	}
	return lo;
}

static void help(void)
{
	printf("Simplistic and fast search of maximal prime gaps 3 to 128 bits\n");
	printf("\n");
	printf("usage : -start xxxxx -end yyyyyy -inc zzzz -radix 10 -h \n");
	printf("\n");
	printf("xxx, yyy and zzz can be decimal numbers or hexadecimal numbers\n");
	printf("like 0x1234567812345678123456780 100000000000000000000  10e20 0b1000010000111\n");
	printf("\n");
	printf("Control of the threads can be done with OMP_NUM_THREADS, and thread affinity with taskset\n");
	printf("\n");
	printf("Examples:\n");
	printf("\n");
	printf("   $ gap_check -start 100 -end 600\n");
	printf("   $ gap_check -start 1e9 -end 1.5e9 -inc 0.2e9\n");
	printf("   $ gap_check -start 0x100000000000000000000 -end 0x1000000000000000000000 -inc 1000000000\n");
	printf("   $ env OMP_NUM_THREADS=4 taskset -c 1,3,5,7 ./tests/gap_check -start 1e20 -end 1e21 -radix 10\n");
	printf("\n");
}

int main(int argc, char **argv)
{
	uint128_t start = (uint128_t) 1 << 68;
	uint128_t end = ((uint128_t) 2 << 68) - 1;
	uint128_t inc = (uint128_t) 1 << 48;
	int radix = 16;
	int digits = 32;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-start")) {
			start = convert128(argv[++i]);
			continue;
		}
		if (!strcmp(argv[i], "-end")) {
			end = convert128(argv[++i]);
			continue;
		}
		if (!strcmp(argv[i], "-inc")) {
			inc = convert128(argv[++i]);
			continue;
		}
		if (!strcmp(argv[i], "-radix")) {
			radix = atoi(argv[++i]);
			if (radix == 16)
				digits = 32;
			else if (radix == 10)
				digits = 40;
			else if (radix == 2)
				digits = 128;
			continue;
		}
		if (!strcmp(argv[i], "-h")) {
			help();
			exit(0);
		}
		help();
		exit(1);
	}

	printf("Simplistic and fast search of maximal prime gaps 3 to 128 bits\n");
	printf("This is test code !\n");
#ifdef __AVX2__
	printf("AVX2 and deterministic primality checks enabled\n");
#else
#error "AVX2 Deterministic primality checks NOT enabled!"
#endif
	printf("\n");
	char tempdisplay[140];
	printf("range from %s\n", display128(tempdisplay, start, radix, digits));
	printf("to         %s\n", display128(tempdisplay, end, radix, digits));
	printf("by steps   %s\n", display128(tempdisplay, inc, radix, digits));
	printf("\n");
	fflush(stdout);

	if (inc > end - start) {
		inc = end - start;
	}

	assert(start >= 7);
	assert(end > start);
	assert(radix == 2 || radix == 10 || radix == 16);

	long t0 = (long)time(NULL);
	uint64_t k = (end - start + inc - 1) / inc;

#pragma omp parallel for
	for (uint64_t x = 0; x < k; x++) {
		uint64_t gap_max = 2;
		uint64_t gap_sz = 0;
		// search between odd primes u and v
		uint128_t u = (start + inc * (uint128_t) x) | 1;	// search lower bound
		uint128_t v = (u + inc) | 1;	// search upper bound
		if (x + 1 != k)
			v += (inc < 10000 ? inc : 10000);	// search over the thread boundary
#ifdef __AVX2__
		_mm256_zeroall();
#endif
		// search the start of the first prime gap
		uint128_t g = search_forwards(u, v);
		uint128_t t = g + 2;
		while (t < v) {
			// search the end of the gap
			t = search_forwards(t, v);
			gap_sz = t - g;
			if (gap_sz >= gap_max) {
				gap_max = gap_sz;
#pragma omp critical
				{
					// this thread has found its own maximal gap,
					// todo : there is no code to centralize that information ,
					// and share to the other threads and processes running the same code.
					// this is only a test program !
					printf("%8ld s.: %s : %5lu\n",
					       (long)time(NULL) - t0, display128(tempdisplay, g, radix, digits), gap_sz);
					fflush(stdout);
				}
			}
			// search the start of the next interesting gap
			g = search_backwards(t, t + gap_max);
			t += gap_max;
		}
	}
}