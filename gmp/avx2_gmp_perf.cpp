#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <x86intrin.h>

#include <gmp.h>
#include "m128_utils.h"
#include "optimized.h"
#include "avx2_sprp.h"

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

	for (unsigned i = 3; i <= 128; i++) {
		uint128_t x = (uint128_t) - 1;
		x /= 7;
		x = ~x;
		x >>= 128 - i;
		x |= 1;

		fflush(stdout);
		fflush(stderr);

		// --------------------------------------------------
		// find a prime number
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
		// check prime numbers
		// --------------------------------------------------
		tavx = 0;
		tgmp = 0;
		for (unsigned j = 0; j < 1000; j++) {
			t = _rdtsc();
			tavx -= t;
			bavx = optimizedSprpTest((uint64_t) x, (uint64_t) (x >> 64)) ? avx2SprpTest((uint64_t) x, (uint64_t) (x >> 64)) : false;
			t = _rdtsc();
			tavx += t;
			tgmp -= t;
			bgmp = (mpz_probab_prime_p(xgmp, 8) != 0);
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
		// average cycles
		// --------------------------------------------------
		double ttavx1 = (double)tavx;
		ttavx1 /= 1000;
		double ttgmp1 = (double)tgmp;
		ttgmp1 /= 1000;

		// --------------------------------------------------
		// find a composite number
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
		// check composite numbers
		// --------------------------------------------------
		tavx = 0;
		tgmp = 0;
		for (unsigned j = 0; j < 1000; j++) {
			t = _rdtsc();
			tavx -= t;
			bavx = optimizedSprpTest((uint64_t) x, (uint64_t) (x >> 64)) ? avx2SprpTest((uint64_t) x, (uint64_t) (x >> 64)) : false;
			t = _rdtsc();
			tavx += t;
			tgmp -= t;
			bgmp = (mpz_probab_prime_p(xgmp, 8) != 0);
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
		// average cycles
		// --------------------------------------------------
		double ttavx2 = (double)tavx;
		ttavx2 /= 1000;
		double ttgmp2 = (double)tgmp;
		ttgmp2 /= 1000;

		// --------------------------------------------------
		// display cycles
		// --------------------------------------------------
		printf("%3d;%9.1f;%9.1f;%9.1f;%9.1f;\n", i, ttavx1, ttgmp1, ttavx2, ttgmp2);
	}
	mpz_clear(xgmp);
	exit(0);
}
