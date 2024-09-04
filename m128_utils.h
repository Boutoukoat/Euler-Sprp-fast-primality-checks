#pragma once

// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// 128 bit utilities
//
// -------------------------------------------------------------------------

#include <stdint.h>
#include <stdbool.h>

typedef unsigned __int128 uint128_t;

// convert a string to 128 bit integer
//
// string format
//
// 0xa1b2c3d4 .... : hexadecimal number 2712847316
// 1234567890 .... : decimal number 1234567890
// 1e12 .......... : decimal number 1000000000000
// 1.5e12 ........ : decimal number 1500000000000
// 0b00010001 .... : binary nnumber 17
// 2^80 .......... : decimal number 1208925819614629174706176
uint128_t convert128(const char *str);

// fill a temp buffer (at least 140 bytes worst case) with 
// 'radix' 2,10,16 representation of 'v', left padded up to 'digits' zeroes.
char *display128(char *temp, uint128_t v, int radix, int digits);

// simple hex display
void my_printf(uint128_t v);
void my_printf(uint64_t v_lo, uint64_t v_hi);

// count leading and trailing zeroes
uint64_t my_clz128(uint128_t v);
uint64_t my_ctz128(uint128_t v);
