// --------------------------------------------------------------------------
//
// check primality tests at random. The reference is gmplib primality test
//
// Stops on failure ..... if ever.
// --------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include <gmp.h>

#include "m128_utils.h"
#include "m_reg.h"
#include "optimized.h"
#include "avx2_sprp.h"

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

	mpz_t xgmp;
	mpz_init(xgmp);

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
			mpz_set_ui(xgmp, (uint64_t) (x >> 64));
			mpz_mul_2exp(xgmp, xgmp, 64);
			mpz_add_ui(xgmp, xgmp, (uint64_t) (x));

			bool bgmp, bopt;

			// reference
			bgmp = mpz_probab_prime_p(xgmp, 8);

			// under test : a sprp test base 2 + enough Sprps tests with small bases (deterministic up to 2^81)
			cycles -= my_rdtsc();
			bopt = optimizedSprpTest(x_lo, x_hi) ? avx2SprpTest(x_lo, x_hi) : false;
			cycles += my_rdtsc();

			check_count += 1;

			if (bgmp != bopt) {
				// some error occured
				// aoutch.  

				printf("Bits ......... : %ld\n", bits);
				printf("Modulus ...... : ");
				my_printf(x);
				printf("\n");
				printf("\n");
				printf("gmp .......... : %s\n", bgmp ? "might be prime" : "composite for sure");
				printf("avx2 ......... : %s\n", bopt ? "might be prime" : "composite for sure");
				printf("\n");
				printf("Fail\n");
				assert(bopt == bgmp);
				abort();
			}
			countPrime += bopt ? 1 : 0;
			countComposite += (!bopt) ? 1 : 0;
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

			printf("%9ld : %16ld primes, %16ld composites %12.1f cycles/primality check (%12.1f nsecs)\n",
			       t2 - t0, countPrime, countComposite, cpf, npf);
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
