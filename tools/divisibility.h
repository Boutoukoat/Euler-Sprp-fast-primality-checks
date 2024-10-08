#pragma once

#include <stdint.h>
#include <m_reg.h>

// ------------------------------------------------------------------------
// fast divisibility of 128 bits numbers by constants 13, 31, 41, 61 .... 
//
// This is probably not the fastest way, modular operations with constant modulus have many tricks (TODO)
// ------------------------------------------------------------------------
bool divisibility_sieve(uint64_t t_lo, uint64_t t_hi);
