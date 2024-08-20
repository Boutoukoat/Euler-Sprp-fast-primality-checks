#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "m128_utils.h"
#include "slow.h"
#include "generic.h"
#include "montgomery.h"
#include "optimized.h"
#include "m_reg.h"

static uint64_t my_random64(void)
{
	// based on linear congruential generator, period = 2^64
	static uint64_t seed = 0x1234567890987654ull;
	seed = seed * 137 + 13;
	// bitwise shuffle
	uint64_t x = seed ^ (seed >> 13) ^ (seed << 13);
	return x;
}

static uint128_t my_random128(void)
{
	uint128_t x = my_random64();
	x ^= 0xaa55aa55aa55aa55ull;
	x <<= 64;
	x += my_random64();
	x ^= 0xc6c6c6c6c6c6c6c6ull;
	return x;
}

int main(int argc __attribute__((unused)), char **argv __attribute__((unused)))
{
	uint64_t start_bit = 3;
	uint64_t end_bit = 128;
	uint64_t delay = 10;	// 10 seconds
	uint64_t t0 = time(NULL);
	uint64_t t1 = time(NULL);
	uint64_t t2 = time(NULL);
	uint64_t seconds = 120;
	uint64_t countPrime = 0;
	uint64_t check_count = 0;
	uint64_t countComposite = 0;
	volatile uint64_t cycles = 0;
	bool fermatTest = true;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-start")) {
			start_bit = strtol(argv[++i], 0, 0);
			continue;
		}
		if (!strcmp(argv[i], "-end")) {
			end_bit = strtol(argv[++i], 0, 0);
			continue;
		}
		if (!strcmp(argv[i], "-time")) {
			seconds = strtol(argv[++i], 0, 0);
			continue;
		}
		if (!strcmp(argv[i], "-sprp")) {
			fermatTest = false;
			continue;
		}
		if (!strcmp(argv[i], "-fermat")) {
			fermatTest = true;
			continue;
		}
		printf("%s -start xxx -end xxxx -time xxxx\n", argv[0]);
		exit(1);
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

	uint128_t start_number = ((uint128_t) 1 << (start_bit - 1)) + 1;
	uint128_t end_number = (end_bit == 128) ? (uint128_t) - 1 : (((uint128_t) 1 << end_bit) - 1);
	printf("Testing numbers in the range [");
	my_printf(start_number);
	printf(", ");
	my_printf(end_number);
	printf("] (bounds included)\n");
	printf("\n");

	fflush(stdout);
	fflush(stderr);

	uint64_t rdtsc_ref_start;
	uint64_t clock_ref_start;
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	clock_ref_start = 1000000000ull * ts.tv_sec + ts.tv_nsec;
	rdtsc_ref_start = my_rdtsc();

	while (t2 - t0 < seconds) {
		for (uint64_t count = 0; count < 1000; count++) {
			// get the bit range at random
			uint64_t bits = my_random64();
			bits = (bits % (end_bit + 1 - start_bit)) + start_bit;	// range 3-128
			uint128_t maskAnd = ((uint128_t) 1 << (bits - 1)) - 1;	// clear bits out of range 
			uint128_t maskOr = ((uint128_t) 1 << (bits - 1)) | 1;	// force msb, force lsb
			//
			// get the number to check
			uint128_t x = my_random128();
			x &= maskAnd;	// clear bits . for 10 bits numbers, mask is 0x1ff
			x |= maskOr;	// force lsb and msb. for 10 bit numbers , mask is 0x201

			uint64_t x_lo, x_hi;
			x_lo = (uint64_t) x;
			x_hi = (uint64_t) (x >> 64);
			bool bgeneric, bopt;

			if (fermatTest) {
				// reference
				bgeneric = genericFermatTest(x_lo, x_hi);

				// under test
				cycles -= my_rdtsc();
				bopt = optimizedFermatTest(x_lo, x_hi);
				cycles += my_rdtsc();
			} else {
				// reference
				bgeneric = montgomerySprpTest(x_lo, x_hi);

				// under test
				cycles -= my_rdtsc();
				bopt = optimizedSprpTest(x_lo, x_hi);
				cycles += my_rdtsc();
			}

			check_count += 1;

			if (bopt != bgeneric) {
				// some error occured
				// aoutch.  

				// get a consensus from different primality algorithms (the error cannot be everywhere)
				bool eslow = my_slowEuler2(x);
				bool f2slow = my_slowFermat2(x);
				bool f3slow = my_slowFermat3(x);
				bool sprp2 = optimizedSprpTest(x_lo, x_hi);
				uint64_t t2 = time(NULL);

				// display the result of many different tests
				printf("Bits ......... : %ld\n", bits);
				printf("Modulus ...... : ");
				my_printf(x);
				printf("\n");
				printf("Fermat-Euler tests\n");
				printf("generic ...... : %s\n", bgeneric ? "might be prime" : "composite for sure");
				printf("optimized .... : %s\n", bopt ? "might be prime" : "composite for sure");
				printf("slow ......... : %s\n", eslow ? "might be prime" : "composite for sure");
				printf("\n");
				printf("Fermat tests\n");
				printf("base 2 ....... : %s\n", f2slow ? "might be prime" : "composite for sure");
				printf("base 3 ....... : %s\n", f3slow ? "might be prime" : "composite for sure");
				printf("\n");
				printf("Sprp tests\n");
				printf("base 2 ....... : %s\n", sprp2 ? "might be prime" : "composite for sure");
				printf("\n");
				printf("%9ld : %16ld pseudoprimes, %16ld composites\n", t2 - t0, countPrime, countComposite);
				printf("Fail\n");
				assert(bopt == bgeneric);
				abort();
			}
			countPrime += bgeneric ? 1 : 0;
			countComposite += (!bgeneric) ? 1 : 0;
		}
		t2 = time(NULL);
		if (t2 - t1 >= delay) {
			uint64_t rdtsc_ref_end;
			uint64_t clock_ref_end;
			clock_gettime(CLOCK_MONOTONIC, &ts);
			clock_ref_end = 1000000000ull * ts.tv_sec + ts.tv_nsec;
			rdtsc_ref_end = my_rdtsc();

			// cycles per fermat
			double cpf = (double)cycles;
			cpf /= countPrime + countComposite;
			// nanosecs per fermat
			double npf = cpf;
			npf *= (double)(clock_ref_end - clock_ref_start);
			npf /= (double)(rdtsc_ref_end - rdtsc_ref_start);

			printf("%9ld : %16ld pseudoprimes, %16ld composites %12.1f cycles/%s (%12.1f nsecs)\n",
			       t2 - t0, countPrime, countComposite, cpf, (fermatTest ? "fermat" : "sprp"), npf);
			t1 = t2;
			fflush(stdout);
			fflush(stderr);
		}
	}
	printf("\n");
	if (check_count) {
		uint64_t rdtsc_ref_end;
		uint64_t clock_ref_end;
		clock_gettime(CLOCK_MONOTONIC, &ts);
		clock_ref_end = 1000000000ull * ts.tv_sec + ts.tv_nsec;	// nsecs
		rdtsc_ref_end = my_rdtsc();	// cycles
		double freq = 1.0;
		freq /= (double)(clock_ref_end - clock_ref_start);	// nsecs
		freq *= 1000000000;	// secs
		freq *= (double)(rdtsc_ref_end - rdtsc_ref_start);	// cycles
		freq /= 1000000;	// Mcycles
		printf("CPU frequency measured %12.1f MHz\n", freq);
	}

	printf("%llu check passed\n", (unsigned long long)check_count);
	fflush(stdout);
	fflush(stderr);
	return 0;
}
