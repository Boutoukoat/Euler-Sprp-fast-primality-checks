
#pragma once
// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// 64-bits utilities
//
// -------------------------------------------------------------------------

#include <stdint.h>

// compute 0 - (1/m) mod (2^64)
// (at run-time, pairs of multiplications are pipelined)
static inline uint64_t montgomeryInverse64(uint64_t mod_lo)
{
	uint64_t x = (3ull * mod_lo) ^ 2ull;	// 5 bits acurate
	uint64_t t = 1ull - mod_lo * x;
	x *= 1 + t;		// 10 bits accurate
	t *= t;
	x *= 1 + t;		// 20 bits accurate
	t *= t;
	x *= 1 + t;		// 40 bits accurate
	t *= t;
	x *= 1 + t;		// 80 bits accurate , i.e. > 64 bits
	return 0 - x;
}

#if 0
// similar code, but slow because all multiplications are sequential
static inline uint64_t slow_montgomeryInverse64(uint64_t mod_lo)
{
uint64_t x = mod_lo + 2 * ((mod_lo +1) & 4) ;
x = x * (2 - mod_lo * x);
x = x * (2 - mod_lo * x);
x = x * (2 - mod_lo * x);
x = x * (2 - mod_lo * x);
return 0 - x;
}
#endif

