#pragma once

// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// low level fast functions
//
// -------------------------------------------------------------------------

#include <stdint.h>

#if defined(__x86_64__)
#include <x86intrin.h>
#endif

typedef unsigned __int128 uint128_t;

// when 0, generate portable C code for any target, including ARM64
// when 1, generate x86_64 instructions
#define INLINE_ASM 1

// hi:lo = a * b
static inline __attribute__((always_inline))
uint64_t my_mul64(uint64_t a, uint64_t b, uint64_t * hi)
{
	if (__builtin_constant_p(hi) && hi == 0) {
		return a * b;
	} else {
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI2__)
		uint64_t lo;
 asm(" mulxq %[a],%[lo],%[hi]\n": [lo] "=r"(lo),[hi] "=r"(*hi): [a] "r"(a), "d"(b):);
		return lo;
#else
		uint64_t lo;
 asm(" mulq %[b]\n": "=a"(lo), "=d"(*hi): [a] "a"(a),[b] "r"(b):"flags");
		return lo;
#endif
#elif (INLINE_ASM && defined(__aarch64__))
		uint64_t lo;

 asm(" mul %[lo], %[a], %[b]\n umulh %[hi], %[a], %[b]\n":[hi] "=&r"(*hi),[lo] "=&r"(lo)
 :		    [a] "r"(a),[b] "r"(b));

		return lo;
#else
		uint128_t tmp = a;
		tmp *= b;
		*hi = (uint64_t) (tmp >> 64);
		return (uint64_t) tmp;
#endif
	}
}

// hi :: lo / d
static inline __attribute__((always_inline))
uint64_t my_div64(uint64_t hi, uint64_t lo, uint64_t d, uint64_t * rem)
{
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(rem) && rem == 0) {
		uint64_t q;
		uint64_t ignore;
 asm(" divq %[d]\n": "=a"(q), "=d"(ignore): "d"(hi), "a"(lo),[d] "r"(d):"flags");
		return q;
	} else {
		uint64_t q;
 asm(" divq %[d]\n": "=a"(q), "=d"(*rem): "d"(hi), "a"(lo),[d] "r"(d):"flags");
		return q;
	}
#else
	uint128_t tmp = ((uint128_t) hi << 64) + lo;
	uint64_t q = (uint64_t) (tmp / d);
	if (__builtin_constant_p(rem) && rem == 0) {
	} else {
		*rem = (uint64_t) (tmp - (uint128_t) q * d);
	}
	return q;
#endif
}

// carry::sum = a + b + carry
static inline __attribute__((always_inline))
uint8_t my_adc64(uint8_t carry_in, uint64_t a, uint64_t b, uint64_t * sum)
{
#if (INLINE_ASM && defined(__x86_64__))
	return _addcarry_u64(carry_in, a, b, (unsigned long long *)sum);
#elif defined(__GNUC__)
	bool c;
	c = __builtin_uaddll_overflow(a, b, (unsigned long long *)sum);
	c |= __builtin_uaddll_overflow(*sum, carry_in, (unsigned long long *)sum);
	return c;
#else
	uint64_t tmp;
	uint8_t carry_out;
	if (__builtin_constant_p(carry_in) && carry_in == 0) {
		if (__builtin_constant_p(a) && a == 0) {
			tmp = b;
			carry_out = 0;
		} else if (__builtin_constant_p(b) && b == 0) {
			tmp = a;
			carry_out = 0;
		} else {
			tmp = a + b;
			carry_out = (tmp < b);
		}
	} else {
		if (__builtin_constant_p(a) && a == 0) {
			tmp = b + carry_in;
			carry_out = (tmp < b);
		} else if (__builtin_constant_p(b) && b == 0) {
			tmp = a + carry_in;
			carry_out = (tmp < a);
		} else {
			tmp = a + carry_in;
			carry_out = (tmp < a);
			tmp += b;
			carry_out |= (tmp < b);
		}
	}
	*sum = tmp;
	return carry_out;
#endif
}

// borrow::diff = a - b - borrow_in
static inline __attribute__((always_inline))
uint8_t my_sbb64(uint8_t borrow_in, uint64_t a, uint64_t b, uint64_t * diff)
{
#if (INLINE_ASM && defined(__x86_64__))
	return _subborrow_u64(borrow_in, a, b, (unsigned long long *)diff);
#elif defined(__GNUC__)
	bool c;
	c = __builtin_usubll_overflow(a, b, (unsigned long long *)diff);
	c |= __builtin_usubll_overflow(*diff, borrow_in, (unsigned long long *)diff);
	return c;
#else
	if (__builtin_constant_p(borrow_in) && borrow_in == 0) {
		if (__builtin_constant_p(a) && a == 0) {
			*diff = -b;
			return 1;
		} else if (__builtin_constant_p(b) && b == 0) {
			*diff = a;
			return 0;
		} else {
			uint64_t tmp = a - b;
			uint8_t borrow = (tmp > a);
			*diff = tmp;
			return borrow;
		}
	} else {
		uint64_t tmp1 = a - borrow_in;
		uint8_t borrow = (tmp1 > a);
		if (__builtin_constant_p(b) && b == 0) {
			*diff = tmp1;
			return borrow;
		} else {
			uint64_t tmp2 = tmp1 - b;
			borrow |= (tmp2 > tmp1);
			*diff = tmp2;
			return borrow;
		}
	}
#endif
}

// (a::b) >>= c
static inline __attribute__((always_inline))
void my_shrd64(uint64_t * a, uint64_t * b, uint64_t c)
{
#if PARANOID
	assert(c < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(c) && c == 0) {
	} else if (__builtin_constant_p(c) && c == 1) {
 asm(" shrdq $1, %1, %0\n shrq $1, %1\n":"=r"(*b), "=r"(*a)
 : "0"(*b), "1"(*a):"flags");
	} else {
 asm(" shrdq %%cl, %1, %0\n shrq %%cl, %1\n":"=&r"(*b), "=&r"(*a)
 : "0"(*b), "1"(*a), "c"(c):"flags");
	}
#else
	if (c == 0) {
	} else {
		(*b) = ((*b) >> c) | ((*a) << (64 - c));
		(*a) >>= c;
	}
#endif
}

// (a::b::c) >>= d
static inline __attribute__((always_inline))
void my_shrd64(uint64_t * a, uint64_t * b, uint64_t * c, uint64_t d)
{
#if PARANOID
	assert(d < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(d) && d == 0) {
	} else if (__builtin_constant_p(d) && d == 1) {
 asm(" shrdq $1, %1, %0\n shrdq $1, %2, %1\n shrq $1, %2\n":"=r"(*c), "=r"(*b), "=r"(*a)
 : "0"(*c), "1"(*b), "2"(*a):"flags");
	} else {
 asm(" shrdq %%cl, %1, %0\n shrdq %%cl, %2, %1\n shrq %%cl, %2\n":"=&r"(*c), "=&r"(*b), "=&r"(*a)
 : "0"(*c), "1"(*b), "2"(*a), "c"(d):"flags");
	}
#else
	if (d == 0) {
	} else {
		(*c) = ((*c) >> d) | ((*b) << (64 - d));
		(*b) = ((*b) >> d) | ((*a) << (64 - d));
		(*a) >>= d;
	}
#endif
}

// (a::b) <<= c
static inline __attribute__((always_inline))
void my_shld64(uint64_t * a, uint64_t * b, uint64_t c)
{
#if PARANOID
	assert(c < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(c) && c == 0) {
	} else if (__builtin_constant_p(c) && c == 1) {
 asm(" addq %0, %0\n adcq %1, %1\n":"=r"(*b), "=r"(*a)
 : "0"(*b), "1"(*a):"flags");
	} else {
 asm(" shldq %%cl, %0, %1\n shlq %%cl, %0\n":"=&r"(*b), "=&r"(*a)
 : "0"(*b), "1"(*a), "c"(c):"flags");
	}
#else
	if (c == 0) {
	} else {
		(*a) = ((*a) << c) | ((*b) >> (64 - c));
		(*b) <<= c;
	}
#endif
}

// (a::b::c) <<= d
static inline __attribute__((always_inline))
void my_shld64(uint64_t * a, uint64_t * b, uint64_t * c, uint64_t d)
{
#if PARANOID
	assert(d < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(d) && d == 0) {
	} else if (__builtin_constant_p(d) && d == 1) {
 asm(" addq %0, %0\n adcq %1, %1\n adcq %2, %2\n":"=r"(*c), "=r"(*b), "=r"(*a)
 : "0"(*c), "1"(*b), "2"(*a):"flags");
	} else {
 asm(" shldq %%cl, %1, %2\n shldq %%cl, %0, %1\n shlq %%cl, %0\n":"=&r"(*c), "=&r"(*b), "=&r"(*a)
 : "0"(*c), "1"(*b), "2"(*a), "c"(d):"flags");
	}
#else
	if (d == 0) {
	} else {
		(*a) = ((*a) << d) | ((*b) >> (64 - d));
		(*b) = ((*b) << d) | ((*c) >> (64 - d));
		(*c) <<= d;
	}
#endif
}

// a >> c
static inline __attribute__((always_inline))
uint64_t my_shr64(uint64_t a, uint64_t c)
{
#if PARANOID
	assert(c < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(c) && c == 0) {
		return a;
	} else {
#ifdef __BMI2__
		uint64_t t;
 asm(" shrxq %2, %1, %0\n": "=r"(t): "r"(a), "r"(c):);
		return t;
#else
 asm(" shrq %%cl, %0\n": "=r"(a): "0"(a), "c"(c):"flags");
		return a;
#endif
	}
#else
	return a >> c;
#endif
}

// a << c
static inline __attribute__((always_inline))
uint64_t my_shl64(uint64_t a, uint64_t c)
{
#if PARANOID
	assert(c < 64);
#endif
#if (INLINE_ASM && defined(__x86_64__))
	if (__builtin_constant_p(c) && c == 0) {
		return a;
	} else if (__builtin_constant_p(c) && c == 1) {
 asm(" addq %0, %0\n": "=r"(a): "0"(a):"flags");
		return a;
	} else {
#ifdef __BMI2__
		uint64_t t;
 asm(" shlxq %2, %1, %0\n": "=r"(t): "r"(a), "r"(c):);
		return t;
#else
 asm(" shlq %%cl, %0\n": "=r"(a): "0"(a), "c"(c):"flags");
		return a;
#endif
	}
#else
	return a << c;
#endif
}

// zero bits after position pos
static inline __attribute__((always_inline))
uint64_t my_bzh64(uint64_t n, uint64_t pos)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI2__
	return _bzhi_u64(n, pos);
#else
	return n & ((1ull << pos) - 1);
#endif
#else
	return n & ((1ull << pos) - 1);
#endif
}

// count trailing zeroes in binary representation 
static inline __attribute__((always_inline))
uint64_t my_ctz32(uint32_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint64_t t;
 asm(" tzcnt %l1, %l0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_ctz(n);
	return 32;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctz(n);
	return 32;
#else
	if (n == 0)
		return 64;
	uint32_t r = 0;
	if ((n & 0xFFFFull) == 0)
		r += 16, n >>= 16;
	if ((n & 0xFFull) == 0)
		r += 8, n >>= 8;
	if ((n & 0xFull) == 0)
		r += 4, n >>= 4;
	if ((n & 0x3ull) == 0)
		r += 2, n >>= 2;
	if ((n & 0x1ull) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count trailing zeroes in binary representation 
static inline __attribute__((always_inline))
uint64_t my_ctz64(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#if defined(__BMI1__)
	uint64_t t;
 asm(" tzcntq %1, %0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_ctzll(n);
	return 64;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_ctzll(n);
	return 64;
#else
	if (n == 0)
		return 64;
	uint64_t r = 0;
	if ((n & 0xFFFFFFFFull) == 0)
		r += 32, n >>= 32;
	if ((n & 0xFFFFull) == 0)
		r += 16, n >>= 16;
	if ((n & 0xFFull) == 0)
		r += 8, n >>= 8;
	if ((n & 0xFull) == 0)
		r += 4, n >>= 4;
	if ((n & 0x3ull) == 0)
		r += 2, n >>= 2;
	if ((n & 0x1ull) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count trailing zeroes in binary representation 
static inline __attribute__((always_inline))
uint64_t my_ctz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_lo) {
		return my_ctz64(n_lo);
	}
	return 64 + my_ctz64(n_hi);
}

// count leading zeroes in binary representation
static inline __attribute__((always_inline))
uint64_t my_clz32(uint32_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint32_t t;
 asm(" lzcnt %l1, %l0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_clz(n);
	return 32;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clz(n);
	return 32;
#else
	if (n == 0)
		return 64;
	uint32_t r = 0;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count leading zeroes in binary representation
static inline __attribute__((always_inline))
uint64_t my_clz64(uint64_t n)
{
#if (INLINE_ASM && defined(__x86_64__))
#ifdef __BMI1__
	uint64_t t;
 asm(" lzcntq %1, %0\n": "=r"(t): "r"(n):"flags");
	return t;
#else
	if (n)
		return __builtin_clzll(n);
	return 64;
#endif
#else
#if defined(__GNUC__)
	if (n)
		return __builtin_clzll(n);
	return 64;
#else
	if (n == 0)
		return 64;
	uint64_t r = 0;
	if ((n & (0xFFFFFFFFull << 32)) == 0)
		r += 32, n <<= 32;
	if ((n & (0xFFFFull << 48)) == 0)
		r += 16, n <<= 16;
	if ((n & (0xFFull << 56)) == 0)
		r += 8, n <<= 8;
	if ((n & (0xFull << 60)) == 0)
		r += 4, n <<= 4;
	if ((n & (0x3ull << 62)) == 0)
		r += 2, n <<= 2;
	if ((n & (0x1ull << 63)) == 0)
		r += 1;
	return r;
#endif
#endif
}

// count leading zeroes in binary representation 
static inline __attribute__((always_inline))
uint64_t my_clz128(uint64_t n_lo, uint64_t n_hi)
{
	if (n_hi) {
		return my_clz64(n_hi);
	}
	return 64 + my_clz64(n_lo);
}

static inline __attribute__((always_inline))
uint64_t my_rdtsc(void)
{
#if defined(__x86_64__)
	// supported by GCC and Clang for x86 platform
	return _rdtsc();
#elif INLINE_ASM && defined(__aarch64__)
	// should be a 64 bits wallclock counter
	// document for old/recent architecture and/or BMC chipsets mention it
	// could be a 56 bit counter.
	uint64_t val;

	asm volatile ("mrs %0, cntvct_el0":"=r" (val));

	// I am not sure what the clock unit is, it depends on pre-scaler setup
	// A multiplication by 32 might be needed on my platform 
	return val * 32;  // aarch64 emulation on x86_64 ?
	return ((val / 3) * 25) << 4;   // maybe for ARM M1 ?
	return val;
#else
#error "todo : unsupported _rdtsc implementation\n"
	return 0;
#endif
}
