
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "m128_utils.h"
#include "montgomery128.cpp"
#include "optimized65.cpp"

uint64_t misc_checks128(uint64_t mod_lo, uint128_t v)
{
	uint64_t check_count = 0;
	uint64_t mod_hi, v_lo, v_hi;
	uint128_t m, e, r;
	mod_hi = 1;
	m = ((uint128_t) mod_hi << 64) + mod_lo;

	// check montgomery inverse
	assert((mod_lo & 1) == 1);
	uint64_t mmagic = montgomeryInverse64(mod_lo);
	check_count += 1;
	if ((mmagic & 1) == 0) {
		printf("modulus ............. : ");
		my_printf(m);
		printf("\n");
		printf("montgomery inverse .. : ");
		my_printf(mmagic);
		printf("\n");
		assert((mmagic & 1) == 1);
		abort();
	}
	// check montgomery inverse inversed
	uint64_t mmagic2 = montgomeryInverse64(mmagic);
	check_count += 1;
	if (mmagic2 != mod_lo) {
		printf("modulus ............. : ");
		my_printf(m);
		printf("\n");
		printf("montgomery inverse .. : ");
		my_printf(mmagic);
		printf("\n");
		printf("montgomery inverse .. : ");
		my_printf(mmagic2);
		printf("\n");
		assert(mmagic2 == mod_lo);
		abort();
	}
	// expected
	v_lo = (uint64_t) v;
	v_hi = (uint64_t) (v >> 64);
	ciosModMul128(&v_lo, &v_hi, v_lo, v_hi, mod_lo, mod_hi, mmagic);
	e = ((uint128_t) v_hi << 64) + v_lo;
	while (e >= m)
		e -= m;

	// test 
	v_lo = (uint64_t) v;
	v_hi = (uint64_t) (v >> 64);
	ciosModSquare(&v_lo, &v_hi, mod_lo, mmagic);
	r = ((uint128_t) v_hi << 64) + v_lo;
	while (r >= m)
		r -= m;

	// result
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
	// expected
	v_lo = (uint64_t) v;
	v_hi = (uint64_t) (v >> 64);
	ciosModMul128(&v_lo, &v_hi, v_lo, v_hi, mod_lo, mod_hi, mmagic);
	ciosModMul128(&v_lo, &v_hi, v_lo, v_hi, mod_lo, mod_hi, mmagic);
	ciosModMul128(&v_lo, &v_hi, v_lo, v_hi, mod_lo, mod_hi, mmagic);
	e = ((uint128_t) v_hi << 64) + v_lo;
	while (e >= m)
		e -= m;

	// test
	v_lo = (uint64_t) v;
	v_hi = (uint64_t) (v >> 64);
	ciosModSquare3(&v_lo, &v_hi, mod_lo, mmagic);
	r = ((uint128_t) v_hi << 64) + v_lo;
	while (r >= m)
		r -= m;

	// result
	check_count += 1;
	if (r != e) {
		printf("modulus ................... : ");
		my_printf(m);
		printf("\n");
		printf("value to square 3 times ... : ");
		my_printf(v);
		printf("\n");
		printf("expected .................. : ");
		my_printf(e);
		printf("\n");
		printf("result .................... : ");
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

	uint128_t v;
	uint64_t mod_lo;

	if (verbose) {
		printf("A few sanity checks\n");
		fflush(stdout);
	}
	// sanity checks
	mod_lo = 0x5c4256770d3e498bull;
	v = (uint128_t) 1ull << 64;
	check_count += misc_checks128(mod_lo, v);

	mod_lo = 0x9d9b054309d65869ull;
	v = 1ull << 63;
	check_count += misc_checks128(mod_lo, v);

	// : 2 smallest modulus
	mod_lo = 1;
	v = ((uint128_t) 1 << 65) + 1;
	check_count += misc_checks128(mod_lo, v);

	mod_lo += 2;
	v = ((uint128_t) 1 << 65) + 3;
	check_count += misc_checks128(mod_lo, v);

	// : 2 largest modulus
	mod_lo = (uint64_t) - 1;
	v = ((uint128_t) 1 << 65) + 1;
	check_count += misc_checks128(mod_lo, v);

	mod_lo -= 2;
	v = ((uint128_t) 1 << 65) + 3;
	check_count += misc_checks128(mod_lo, v);

	if (verbose) {
		printf("Multiple number combination\n");
		fflush(stdout);
	}
	// check multiple combinations for m < 65 bits
	for (mod_lo = 3; mod_lo <= 0xffffffffffffffffull / 3; mod_lo *= 3) {
		// mod_lo is always odd
		for (v = 0; v <= mod_lo; v = v * 7 + 1) {
			check_count += misc_checks128(mod_lo, v);
		}
	}

	assert(check_count != 0);
	printf("%lu checks passed\n", check_count);
	return 0;
}
