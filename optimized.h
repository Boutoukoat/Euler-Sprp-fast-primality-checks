
#pragma once
// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// fastest code for primality testing 
//
// -------------------------------------------------------------------------

#include <stdint.h>
#include <stdbool.h>

// fastest code for 65 bits operations
bool optimizedSprpTest65(uint64_t n_lo);
bool optimizedFermatTest65(uint64_t n_lo);

// preferred entry points for all values of n_hi and n_lo
// (this should be the fastest code)
bool optimizedSprpTest(uint64_t n_lo, uint64_t n_hi);
bool optimizedFermatTest(uint64_t n_lo, uint64_t n_hi);
