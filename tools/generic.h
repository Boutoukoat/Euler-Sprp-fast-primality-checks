#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef unsigned __int128 uint128_t;

bool genericFermatTest(uint64_t n_lo, uint64_t n_hi);
bool genericFermatTest(uint64_t base, uint64_t n_lo, uint64_t n_hi);
