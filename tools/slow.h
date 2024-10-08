#pragma once

#include <stdint.h>
#include <stdbool.h>

typedef unsigned __int128 uint128_t;

// shift and add 128 bits
uint128_t my_slowReduce(uint128_t x, uint128_t mod);

// based on shift and add 128 bits
uint128_t my_slowModSqr(uint128_t x, uint128_t mod);

// based on shift and add 128 bits
bool my_slowSprp2(uint128_t mod);
bool my_slowSprp3(uint128_t mod);
bool my_slowSprp5(uint128_t mod);
bool my_slowFermat2(uint128_t mod);
bool my_slowFermat3(uint128_t mod);
bool my_slowFermat5(uint128_t mod);
bool my_slowEuler2(uint128_t mod);
