#pragma once

#include "stdint.h"
#include "stdbool.h"

typedef unsigned __int128 uint128_t;

// ----------------------------------------------------------------
//
// Internal entry points
//
// v : number to check > 2 bits and less than 2^128
// mm : -1/v mod 2^32   (montgomery inverse)
// on : 2^(size_v) mod v (montgomery one)  size_v is 32,64,96 or128
// bases : exactly 4 bases (montgomery domain)
//         test 4 sprps bases 11,13,17,19 means
//         bases[0] = (11 * one)) % v
//         bases[1] = (13 * one) % v
//         bases[2] = (17 * one) % v
//         bases[3] = (19 * one) % v
//
// return true : might be prime (depends on set of bases, and on v)
// return false : composite for sure
//
// Choice of bases is to the user
//
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
// https://miller-rabin.appspot.com/
//
// ----------------------------------------------------------------
bool avx2_sprp1(uint32_t v, uint32_t mm, uint32_t on, uint32_t * bases);	// 2-31 bits
bool avx2_sprp2(uint64_t v, uint32_t mm, uint64_t on, uint64_t * bases);	// 32-63 bits
bool avx2_sprp3(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases);	// 64-95 bits
bool avx2_sprp4(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases);	// 96-127 bits

// ----------------------------------------------------------------
//
// avx2-based deterministic primality check
//
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Testing_against_small_sets_of_bases
//
// v : number to check > 2 bits and less than 2^128
//
// return true : prime for sure   (up to now, proven under 81.5 bits)
// return false : composite for sure
//
// ----------------------------------------------------------------
bool avx2SprpTest(uint64_t v_lo, uint64_t v_hi);
