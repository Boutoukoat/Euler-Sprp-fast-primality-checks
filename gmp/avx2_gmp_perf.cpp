// -------------------------------------------------------------------------------------
//
// Fast primality tests : compare GMP and multiple sprp tests
// 
// Deterministic to 81.5 bits https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
// and better choices of bases could be taken from https://miller-rabin.appspot.com/
// -------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>

#include <gmp.h>
#include "m128_utils.h"
#include "m_reg.h"
#include "optimized.h"
#include "avx2_sprp.h"

// little sieve for small constants and large numbers
#include "divisibility.h"

int main(int argc, char **argv)
{
	bool verbose = false;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-v")) {
			verbose = true;
			continue;
		}
		printf("Wrong argument on command line : %s\n", argv[i]);
		exit(1);
	}

	bool bavx, bgmp;
	uint64_t t;
	volatile uint64_t tavx, tgmp;
	mpz_t xgmp;
	mpz_init(xgmp);

	// CSV header
	printf("\"Bits\";");
	printf("\"Prime AVX2\";");
	printf("\"Prime GMP\";");
	printf("\"Composite AVX2\";");
	printf("\"Composite GMP\";");
	printf("\n");

	for (unsigned i = 3; i <= 128; i++) {
		// fill a mantissa with 0 and 1
		uint128_t x = (uint128_t) - 1;
		x /= 7;
		// make sure a few msbbits are set
		x = ~x;
		// adjust the size of the number
		x >>= 128 - i;
		// make it odd
		x |= 1;

		fflush(stdout);
		fflush(stderr);

		// --------------------------------------------------
		// search a prime number, brute force
		// --------------------------------------------------
		do {
			x += 2;
			mpz_set_ui(xgmp, (uint64_t) (x >> 64));
			mpz_mul_2exp(xgmp, xgmp, 64);
			mpz_add_ui(xgmp, xgmp, (uint64_t) (x));
		}
		while (!mpz_probab_prime_p(xgmp, 16));

		if (verbose) {
			printf("Testing prime ");
			my_printf(x);
			printf("\n");
		}
		// --------------------------------------------------
		// check the prime number
		// --------------------------------------------------
		tavx = 0;
		tgmp = 0;
		for (unsigned j = 0; j < 1000; j++) {
			t = _rdtsc();
			tavx -= t;
			bavx = (divisibility_sieve((uint64_t) x, (uint64_t) (x >> 64))
				&& optimizedSprpTest((uint64_t) x, (uint64_t) (x >> 64))
				&& avxSprpTest((uint64_t) x, (uint64_t) (x >> 64)));
			t = _rdtsc();
			tavx += t;
			tgmp -= t;
			// between 3 and 20 rounds
			bgmp = !!mpz_probab_prime_p(xgmp, 3 + (17 * i) / 128);
			t = _rdtsc();
			tgmp += t;
			if (bavx != bgmp) {
				printf("Error at prime : ");
				my_printf(x);
				printf("\n");
				exit(1);
			}
		}

		// --------------------------------------------------
		// get average cycles
		// --------------------------------------------------
		double ttavx1 = (double)tavx;
		ttavx1 /= 1000;
		double ttgmp1 = (double)tgmp;
		ttgmp1 /= 1000;

		// --------------------------------------------------
		// search a composite number, brute force
		// --------------------------------------------------
		do {
			x += 2;
			mpz_set_ui(xgmp, (uint64_t) (x >> 64));
			mpz_mul_2exp(xgmp, xgmp, 64);
			mpz_add_ui(xgmp, xgmp, (uint64_t) (x));
		}
		while (mpz_probab_prime_p(xgmp, 16));

		if (verbose) {
			printf("Testing composite ");
			my_printf(x);
			printf("\n");
		}
		// --------------------------------------------------
		// check the composite number
		// --------------------------------------------------
		tavx = 0;
		tgmp = 0;
		for (unsigned j = 0; j < 1000; j++) {
			t = _rdtsc();
			tavx -= t;
			bavx = (divisibility_sieve((uint64_t) x, (uint64_t) (x >> 64))
				&& optimizedSprpTest((uint64_t) x, (uint64_t) (x >> 64))
				&& avxSprpTest((uint64_t) x, (uint64_t) (x >> 64)));
			t = _rdtsc();
			tavx += t;
			tgmp -= t;
			// between 3 and 20 rounds
			bgmp = !!mpz_probab_prime_p(xgmp, 3 + (17 * i) / 128);
			t = _rdtsc();
			tgmp += t;
			if (bavx != bgmp) {
				printf("Error at composite : ");
				my_printf(x);
				printf("\n");
				exit(1);
			}
		}

		// --------------------------------------------------
		// get average cycles
		// --------------------------------------------------
		double ttavx2 = (double)tavx;
		ttavx2 /= 1000;
		double ttgmp2 = (double)tgmp;
		ttgmp2 /= 1000;

		// --------------------------------------------------
		// display cycles  (CSV-compatible format)
		// --------------------------------------------------
		printf("%3d;%9.1f;%9.1f;%9.1f;%9.1f;\n", i, ttavx1, ttgmp1, ttavx2, ttgmp2);
	}
	mpz_clear(xgmp);
	exit(0);
}
