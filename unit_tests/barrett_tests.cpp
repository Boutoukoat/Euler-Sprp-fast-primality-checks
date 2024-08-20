
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "m128_utils.h"
#include "slow.h"
#include "barrett.cpp"

static uint128_t barrett_square_reduce128(uint128_t v, uint64_t mod_lo, uint64_t mod_hi)
{
	uint64_t barrett_lo, barrett_hi;

	reciprocal(mod_lo, mod_hi, &barrett_lo, &barrett_hi);
	return barrett_mod_square(v, mod_lo, mod_hi, barrett_lo, barrett_hi);
}

static uint64_t barrett_square_reduce64(uint64_t v, uint64_t mod_lo)
{
	uint64_t barrett_lo, barrett_hi;

	reciprocal64(mod_lo, &barrett_lo, &barrett_hi);
	return barrett_mod_square64(v, mod_lo, barrett_lo, barrett_hi);
}

uint64_t misc_checks64(uint64_t m, uint64_t v)
{
	uint64_t check_count = 0;
	uint128_t e, r;

	// expected
	e = my_slowModSqr(v, m);
	// result
	r = barrett_square_reduce64(v, m);
	check_count += 1;
	if (r != e) {
		printf("modulus ............. : ");
		my_printf(m);
		printf("\n");
		printf("value to square ..... : ");
		my_printf(v);
		printf("\n");
		printf("expected ............ : ");
		my_printf(e);
		printf("\n");
		printf("result .............. : ");
		my_printf(r);
		printf("\n");
		assert(r == e);
		abort();
	}
	return check_count;

}

uint64_t misc_checks128(uint128_t m, uint128_t v)
{
	uint64_t check_count = 0;
	uint64_t mod_hi, mod_lo;
	uint128_t e, r;
	mod_lo = (uint64_t) m;
	mod_hi = (uint64_t) (m >> 64);

	// expected
	e = my_slowModSqr(v, m);
	// result
	r = barrett_square_reduce128(v, mod_lo, mod_hi);
	check_count += 1;
	if (r != e) {
		printf("modulus ............. : ");
		my_printf(m);
		printf("\n");
		printf("value to square ..... : ");
		my_printf(v);
		printf("\n");
		printf("expected ............ : ");
		my_printf(e);
		printf("\n");
		printf("result .............. : ");
		my_printf(r);
		printf("\n");
		assert(r == e);
		abort();
	}
	return check_count;

}

int main(int argc, char **argv)
{
	uint64_t check_count = 0;
	bool verbose = false;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-v") || !strcmp(argv[i], "-verbose")) {
			verbose = true;
		}
	}

	uint128_t m, v;
	uint64_t mod_lo, mod_hi, u;

	if (verbose) {
		printf("A few sanity checks\n");
		fflush(stdout);
	}
	// 64 bits 
	mod_lo = 0x5c4256770d3e498bull;
	u = 1ull << 53;
	check_count += misc_checks64(mod_lo, u);

	// 64 bits 
	mod_lo = 0x9d9b054309d65868ull;
	u = 1ull << 53;
	check_count += misc_checks64(mod_lo, u);

	// 
	m = ((uint128_t) 0x9d9b054309d65868ull << 64) + 0x5c4256770d3e498bull;
	v = (uint128_t) 1ull << 78;
	check_count += misc_checks128(m, v);

	// 
	m = (uint128_t) 0x33333333333333ull *0x55555555555555ull;
	v = (uint128_t) 0x11111111111111ull *0x77777777777777ull;
	check_count += misc_checks128(m, v);

	// : 2 smallest modulus
	m = ((uint128_t) 1 << 64) + 1;
	v = m * 2 + 1;
	check_count += misc_checks128(m, v);

	m += 2;
	v = m * 2 + 1;
	check_count += misc_checks128(m, v);

	// : 2 largest modulus
	m = (uint128_t) - 1;
	v = m - 1;
	check_count += misc_checks128(m, v);

	m -= 2;
	v = m + 1;
	check_count += misc_checks128(m, v);

	// : smallest 128 bit modulus
	m = ((uint128_t) 1 << 127) + 1;
	v = m + 100;
	check_count += misc_checks128(m, v);

	if (verbose) {
		printf("65-128 bits checks\n");
		fflush(stdout);
	}
	// check multiple combinations for m > 64 bits
	for (mod_hi = 1; mod_hi <= (0xffffffffffffffffull - 1) / 3; mod_hi = mod_hi * 3 + 1) {
		// mod_hi is > 0 and alternates odd/even
		for (mod_lo = 1; mod_lo <= 0xffffffffffffffffull / 5; mod_lo *= 5) {
			// mod_lo is always odd
			m = ((uint128_t) mod_hi << 64) + mod_lo;
			for (v = 0; v <= (m - 1) / 7; v = v * 7 + 1) {
				check_count += misc_checks128(m, v);
			}
		}
	}

	if (verbose) {
		printf("3-64 bits checks\n");
		fflush(stdout);
	}
	// check multiple combinations for m < 65 bits
	for (mod_lo = 5; mod_lo <= 0xffffffffffffffffull / 5; mod_lo *= 5) {
		// mod_lo is always odd
		for (u = 0; u <= mod_lo / 5; u = u * 3 + 1) {
			check_count += misc_checks64(mod_lo, u);
		}
	}

	assert(check_count != 0);
	printf("%lu checks passed\n", check_count);
	return 0;
}
