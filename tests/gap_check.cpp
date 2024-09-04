#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#if PARANOID
#include "montgomery.h"
#endif

#include "optimized.h"

// utilities to convert strings from/to 128 bits inumbers
#include "m128_utils.h"

// utilities for divisibility sieve
#include "m_reg.h"

// use a simplistic wheel sieve mod 2,3,5,7,11
#include "gap_check.wheel"

// simple divisibility tests for large 128-bit numbers and larger
#include "divisibility.h"

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
		bool bs = divisibility_sieve(t_lo, t_hi);
		if (bs) {
#if PARANOID
			bool bg = montgomerySprpTest(t_lo, t_hi);
			bool ba = optimizedSprpTest(t_lo, t_hi);
			if (ba != bg) {
				printf("0x%16.16lx%16.16lx\n", t_hi, t_lo);
				assert(ba == bg);
			}
#else
			bool ba = optimizedSprpTest(t_lo, t_hi);
#endif
			if (ba) {
				// pseudoprime found
				break;
			}
		}
		// move backwards to the next sieve candidate
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
		bool bs = divisibility_sieve(t_lo, t_hi);
		if (bs) {
#if PARANOID
			bool bg = montgomerySprpTest(t_lo, t_hi);
			bool ba = optimizedSprpTest(t_lo, t_hi);
			if (ba != bg) {
				printf("0x%16.16lx%16.16lx\n", t_hi, t_lo);
				assert(ba == bg);
			}
#else
			bool ba = optimizedSprpTest(t_lo, t_hi);
#endif
			if (ba) {
				// pseudoprime found
				break;
			}
		}
		// move to the next sieve candidate
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
					// this thread has found its own maximal gap
					// todo : there is no code to centralize that information 
					// and share to the other threads and processes running the same code
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
