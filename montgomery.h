
#pragma once
// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// Entry points for primality tests based on montgomery reduction
//
// -------------------------------------------------------------------------

#include <stdint.h>
#include <stdbool.h>

// internal entry points
bool montgomeryFermatTest128(uint64_t n_lo, uint64_t n_hi);
bool montgomeryFermatTest64(uint64_t n_lo);
bool montgomeryFermatTest8(uint64_t n_lo);
bool montgomerySprpTest128(uint64_t n_lo, uint64_t n_hi);
bool montgomerySprpTest64(uint64_t n_lo);

// preferred entry points for all values of n_hi and n_lo
// (this is the textbook variant, not the fastest)
bool montgomeryFermatTest(uint64_t n_lo, uint64_t n_hi);
bool montgomerySprpTest(uint64_t n_lo, uint64_t n_hi);
