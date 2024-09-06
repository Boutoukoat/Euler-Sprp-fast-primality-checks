
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <x86intrin.h>

#include "avx2_sprp.h"
#include "m_reg.h"
#include "m128_utils.h"
#include "m64_utils.h"

#ifdef __AVX2__

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract1(__m256i * mr, __m256i mp[1])
{
	__m256i t[1], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract21(__m256i * mr, __m256i mc, __m256i mp[1])
{
	__m256i t[1], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mc);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_subtract1(__m256i * mr, __m256i ma[1], __m256i mb[1])
{
	mr[0] = _mm256_sub_epi64(ma[0], mb[0]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_cmp_next1(__m256i * mask, __m256i ma[1], __m256i mb[1])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_neg_cmp_next1(__m256i * mask, __m256i ma[1], __m256i mb[1])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modmul1(__m256i * mr, __m256i mb[1], __m256i mp[1], __m256i mmagic)
{
	__m256i t[3];
	__m256i m, c, cs, d, ds, a;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #2
	t[1] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[0]);	// #3
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, t[1]);
	mr[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	ds = _mm256_srli_epi64(ds, 32);
	avx2_modsubtract21(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsqu1(__m256i * mr, __m256i mp[1], __m256i mmagic)
{
	__m256i t[3];
	__m256i a, m, c, c1, c2, cs, ct, cd;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, a);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #2
	t[1] = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[0]);	// #3
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[1]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_srli_epi64(cs, 32);
	avx2_modsubtract21(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx2_sprp1(uint32_t v, uint32_t mm, uint32_t on, uint32_t * bases)
{
	__m256i p[1], r[1], one[1], m1[1], b[1];
	uint64_t bit, k;
	uint32_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 0),
				0, (uint32_t) (bases[1] >> 0), 0, (uint32_t) (bases[2] >> 0), 0, (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
// p - 1
	avx2_subtract1(m1, p, one);
// first value
	r[0] = b[0];
// MR exponentiation bit per bit
	k = my_ctz32(v - 1);
	s = v >> k;
	bit = 31 - my_clz32(s);
	while (bit > 0) {
		bit--;
		// Square
		avx2_modsqu1(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx2_modmul1(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[1], mask = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	if (avx2_cmp_next1(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx2_cmp_next1(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx2_cmp_next1(&mask, r, m1)) {
			return true;
		}
		// square
		avx2_modsqu1(r, p, mmagic);
		// check current result == 1
		if (avx2_neg_cmp_next1(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx2_cmp_next1(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract2(__m256i * mr, __m256i mp[2])
{
	__m256i t[2], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract32(__m256i * mr, __m256i mc, __m256i mp[2])
{
	__m256i t[2], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mc);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_subtract2(__m256i * mr, __m256i ma[2], __m256i mb[2])
{
	__m256i cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[1]);
	mr[1] = _mm256_sub_epi64(cs, mb[1]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_cmp_next2(__m256i * mask, __m256i ma[2], __m256i mb[2])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_neg_cmp_next2(__m256i * mask, __m256i ma[2], __m256i mb[2])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modmul2(__m256i * mr, __m256i mb[2], __m256i mp[2], __m256i mmagic)
{
	__m256i t[4];
	__m256i m, c, cs, d, ds, a;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm256_mul_epu32(m, mp[0]);	// #3
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #4
	cs = _mm256_add_epi64(cs, c);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#5
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, c);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	t[2] = _mm256_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm256_mul_epu32(a, mb[0]);	// #6
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #7
	ds = _mm256_mul_epu32(m, mp[0]);	// #8
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #9
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#10
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	mr[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[3] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[2]);
	mr[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, t[3]);

	avx2_modsubtract32(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsqu2(__m256i * mr, __m256i mp[2], __m256i mmagic)
{
	__m256i t[4];
	__m256i a, m, c, c1, c2, cs, ct, cd;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, a);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[1]);	// #2
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[3] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #3
	cs = _mm256_mul_epu32(m, mp[0]);	// #4
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #5
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[2] = _mm256_add_epi64(c, t[3]);

	a = mr[1];
	cs = _mm256_mul_epu32(a, a);	// #6
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c1, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[3] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #7
	cs = _mm256_mul_epu32(m, mp[0]);	// #8
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #9
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[2]);
	mr[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[3]);

	avx2_modsubtract32(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx2_sprp2(uint64_t v, uint32_t mm, uint64_t on, uint64_t * bases)
{
	__m256i p[2], r[2], one[2], m1[2], b[2];
	uint64_t bit, k;
	uint64_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 0),
				0, (uint32_t) (bases[1] >> 0), 0, (uint32_t) (bases[2] >> 0), 0, (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 32),
				0, (uint32_t) (bases[1] >> 32), 0, (uint32_t) (bases[2] >> 32), 0, (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
// p - 1
	avx2_subtract2(m1, p, one);
// first value
	r[0] = b[0];
	r[1] = b[1];
// MR exponentiation bit per bit
	k = my_ctz64(v - 1);
	s = v >> k;
	bit = 63 - my_clz64(s);
	while (bit > 0) {
		bit--;
		// Square
		avx2_modsqu2(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx2_modmul2(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[2], mask = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	if (avx2_cmp_next2(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx2_cmp_next2(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx2_cmp_next2(&mask, r, m1)) {
			return true;
		}
		// square
		avx2_modsqu2(r, p, mmagic);
		// check current result == 1
		if (avx2_neg_cmp_next2(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx2_cmp_next2(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract3(__m256i * mr, __m256i mp[3])
{
	__m256i t[3], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[2]);
	cs = _mm256_sub_epi64(cs, mp[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
	mr[2] = _mm256_blendv_epi8(t[2], mr[2], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract43(__m256i * mr, __m256i mc, __m256i mp[3])
{
	__m256i t[3], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[2]);
	cs = _mm256_sub_epi64(cs, mp[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mc);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
	mr[2] = _mm256_blendv_epi8(t[2], mr[2], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_subtract3(__m256i * mr, __m256i ma[3], __m256i mb[3])
{
	__m256i cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[1]);
	cs = _mm256_sub_epi64(cs, mb[1]);
	mr[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[2]);
	mr[2] = _mm256_sub_epi64(cs, mb[2]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_cmp_next3(__m256i * mask, __m256i ma[3], __m256i mb[3])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_neg_cmp_next3(__m256i * mask, __m256i ma[3], __m256i mb[3])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modmul3(__m256i * mr, __m256i mb[3], __m256i mp[3], __m256i mmagic)
{
	__m256i t[5];
	__m256i m, c, cs, d, ds, a;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm256_mul_epu32(m, mp[0]);	// #3
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #4
	cs = _mm256_add_epi64(cs, c);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#5
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #6
	cs = _mm256_add_epi64(cs, c);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#7
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, c);
	t[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	t[3] = _mm256_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm256_mul_epu32(a, mb[0]);	// #8
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #9
	ds = _mm256_mul_epu32(m, mp[0]);	// #10
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #11
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#12
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #13
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#14
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[4] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[3]);
	t[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	t[3] = _mm256_add_epi64(d, t[4]);

	a = mr[2];
	cs = _mm256_mul_epu32(a, mb[0]);	// #15
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #16
	ds = _mm256_mul_epu32(m, mp[0]);	// #17
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #18
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#19
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	mr[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #20
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#21
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	mr[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[4] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[3]);
	mr[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, t[4]);

	avx2_modsubtract43(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsqu3(__m256i * mr, __m256i mp[3], __m256i mmagic)
{
	__m256i t[5];
	__m256i a, m, c, c1, c2, cs, ct, cd;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, a);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[1]);	// #2
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[2]);	// #3
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	cs = _mm256_add_epi64(cs, c2);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[4] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #4
	cs = _mm256_mul_epu32(m, mp[0]);	// #5
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #6
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #7
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[3]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[3] = _mm256_add_epi64(c, t[4]);

	a = mr[1];
	cs = _mm256_mul_epu32(a, a);	// #8
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[2]);	// #9
	cs = _mm256_add_epi64(ct, c1);
	cs = _mm256_add_epi64(cs, t[2]);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	cs = _mm256_add_epi64(cs, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[4] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #10
	cs = _mm256_mul_epu32(m, mp[0]);	// #11
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #12
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #13
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[3]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[3] = _mm256_add_epi64(c, t[4]);

	a = mr[2];
	cs = _mm256_mul_epu32(a, a);	// #14
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c1, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[4] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #15
	cs = _mm256_mul_epu32(m, mp[0]);	// #16
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #17
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #18
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	mr[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[3]);
	mr[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[4]);

	avx2_modsubtract43(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx2_sprp3(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m256i p[3], r[3], one[3], m1[3], b[3];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 0),
				0, (uint32_t) (bases[1] >> 0), 0, (uint32_t) (bases[2] >> 0), 0, (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 32),
				0, (uint32_t) (bases[1] >> 32), 0, (uint32_t) (bases[2] >> 32), 0, (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b[2] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 64),
				0, (uint32_t) (bases[1] >> 64), 0, (uint32_t) (bases[2] >> 64), 0, (uint32_t) (bases[3] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
// p - 1
	avx2_subtract3(m1, p, one);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
// MR exponentiation bit per bit
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 0) {
		bit--;
		// Square
		avx2_modsqu3(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx2_modmul3(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[3], mask = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	zero[2] = _mm256_setzero_si256();
	if (avx2_cmp_next3(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx2_cmp_next3(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx2_cmp_next3(&mask, r, m1)) {
			return true;
		}
		// square
		avx2_modsqu3(r, p, mmagic);
		// check current result == 1
		if (avx2_neg_cmp_next3(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx2_cmp_next3(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract4(__m256i * mr, __m256i mp[4])
{
	__m256i t[4], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[2]);
	cs = _mm256_sub_epi64(cs, mp[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[3]);
	cs = _mm256_sub_epi64(cs, mp[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
	mr[2] = _mm256_blendv_epi8(t[2], mr[2], cs);
	mr[3] = _mm256_blendv_epi8(t[3], mr[3], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsubtract54(__m256i * mr, __m256i mc, __m256i mp[4])
{
	__m256i t[4], cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(mr[0], mp[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[1]);
	cs = _mm256_sub_epi64(cs, mp[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[2]);
	cs = _mm256_sub_epi64(cs, mp[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mr[3]);
	cs = _mm256_sub_epi64(cs, mp[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, mc);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	mr[0] = _mm256_blendv_epi8(t[0], mr[0], cs);
	mr[1] = _mm256_blendv_epi8(t[1], mr[1], cs);
	mr[2] = _mm256_blendv_epi8(t[2], mr[2], cs);
	mr[3] = _mm256_blendv_epi8(t[3], mr[3], cs);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_subtract4(__m256i * mr, __m256i ma[4], __m256i mb[4])
{
	__m256i cs, zero = _mm256_setzero_si256();
	cs = _mm256_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[1]);
	cs = _mm256_sub_epi64(cs, mb[1]);
	mr[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[2]);
	cs = _mm256_sub_epi64(cs, mb[2]);
	mr[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	cs = _mm256_shuffle_epi32(cs, _MM_SHUFFLE(3, 3, 1, 1));
	cs = _mm256_add_epi64(cs, ma[3]);
	mr[3] = _mm256_sub_epi64(cs, mb[3]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_cmp_next4(__m256i * mask, __m256i ma[4], __m256i mb[4])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[3], mb[3]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_neg_cmp_next4(__m256i * mask, __m256i ma[4], __m256i mb[4])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi32(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi32(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi32(ma[3], mb[3]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modmul4(__m256i * mr, __m256i mb[4], __m256i mp[4], __m256i mmagic)
{
	__m256i t[6];
	__m256i m, c, cs, d, ds, a;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm256_mul_epu32(m, mp[0]);	// #3
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #4
	cs = _mm256_add_epi64(cs, c);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#5
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #6
	cs = _mm256_add_epi64(cs, c);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#7
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[3]);	// #8
	cs = _mm256_add_epi64(cs, c);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[3]);	//#9
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[3]);
	t[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, c);
	t[3] = _mm256_blend_epi32(ds, zero, 0xaa);
	t[4] = _mm256_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm256_mul_epu32(a, mb[0]);	// #10
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #11
	ds = _mm256_mul_epu32(m, mp[0]);	// #12
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #13
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#14
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #15
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#16
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[3]);	// #17
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[3]);	//#18
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[3]);
	t[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[4]);
	t[3] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	t[4] = _mm256_add_epi64(d, t[5]);

	a = mr[2];
	cs = _mm256_mul_epu32(a, mb[0]);	// #19
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #20
	ds = _mm256_mul_epu32(m, mp[0]);	// #21
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #22
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#23
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	t[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #24
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#25
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	t[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[3]);	// #26
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[3]);	//#27
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[3]);
	t[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[4]);
	t[3] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	t[4] = _mm256_add_epi64(d, t[5]);

	a = mr[3];
	cs = _mm256_mul_epu32(a, mb[0]);	// #28
	cs = _mm256_add_epi64(cs, t[0]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	m = _mm256_mul_epu32(mmagic, t[0]);	// #29
	ds = _mm256_mul_epu32(m, mp[0]);	// #30
	ds = _mm256_add_epi64(ds, t[0]);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[1]);	// #31
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[1]);	//#32
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[1]);
	mr[0] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[2]);	// #33
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[2]);	//#34
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[2]);
	mr[1] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_mul_epu32(a, mb[3]);	// #35
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	ds = _mm256_mul_epu32(m, mp[3]);	//#36
	ds = _mm256_add_epi64(ds, d);
	ds = _mm256_add_epi64(ds, t[3]);
	mr[2] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);
	ds = _mm256_add_epi64(d, t[4]);
	mr[3] = _mm256_blend_epi32(ds, zero, 0xaa);
	d = _mm256_srli_epi64(ds, 32);
	ds = _mm256_add_epi64(d, t[5]);

	avx2_modsubtract54(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx2_modsqu4(__m256i * mr, __m256i mp[4], __m256i mmagic)
{
	__m256i t[6];
	__m256i a, m, c, c1, c2, cs, ct, cd;
	__m256i zero = _mm256_setzero_si256();

	a = mr[0];
	cs = _mm256_mul_epu32(a, a);	// #1
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[1]);	// #2
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[2]);	// #3
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	cs = _mm256_add_epi64(cs, c2);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[3]);	// #4
	cs = _mm256_add_epi64(ct, c1);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	cs = _mm256_add_epi64(cs, c2);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #5
	cs = _mm256_mul_epu32(m, mp[0]);	// #6
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #7
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #8
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[3]);	// #9
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[4] = _mm256_add_epi64(c, t[5]);

	a = mr[1];
	cs = _mm256_mul_epu32(a, a);	// #10
	cs = _mm256_add_epi64(cs, t[1]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[2]);	// #11
	cs = _mm256_add_epi64(ct, c1);
	cs = _mm256_add_epi64(cs, t[2]);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[3]);	// #12
	cs = _mm256_add_epi64(ct, c1);
	cs = _mm256_add_epi64(cs, t[3]);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	cs = _mm256_add_epi64(cs, c2);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	cs = _mm256_add_epi64(cs, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #13
	cs = _mm256_mul_epu32(m, mp[0]);	// #14
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #15
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #16
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[3]);	// #17
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[4] = _mm256_add_epi64(c, t[5]);

	a = mr[2];
	cs = _mm256_mul_epu32(a, a);	// #18
	cs = _mm256_add_epi64(cs, t[2]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	ct = _mm256_mul_epu32(a, mr[3]);	// #19
	cs = _mm256_add_epi64(ct, c1);
	cs = _mm256_add_epi64(cs, t[3]);
	cd = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(ct, cd);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c2 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c2, c1);
	cs = _mm256_add_epi64(cs, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #20
	cs = _mm256_mul_epu32(m, mp[0]);	// #21
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #22
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	t[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #23
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	t[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[3]);	// #24
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	t[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[4]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	t[4] = _mm256_add_epi64(c, t[5]);

	a = mr[3];
	cs = _mm256_mul_epu32(a, a);	// #25
	cs = _mm256_add_epi64(cs, t[3]);
	t[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c1, t[4]);
	t[4] = _mm256_blend_epi32(cs, zero, 0xaa);
	t[5] = _mm256_srli_epi64(cs, 32);

	m = _mm256_mul_epu32(mmagic, t[0]);	// #26
	cs = _mm256_mul_epu32(m, mp[0]);	// #27
	cs = _mm256_add_epi64(cs, t[0]);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[1]);	// #28
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[1]);
	mr[0] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[2]);	// #29
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[2]);
	mr[1] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_mul_epu32(m, mp[3]);	// #30
	cs = _mm256_add_epi64(cs, c);
	cs = _mm256_add_epi64(cs, t[3]);
	mr[2] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[4]);
	mr[3] = _mm256_blend_epi32(cs, zero, 0xaa);
	c = _mm256_srli_epi64(cs, 32);
	cs = _mm256_add_epi64(c, t[5]);

	avx2_modsubtract54(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx2_sprp4(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m256i p[4], r[4], one[4], m1[4], b[4];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 0),
				0, (uint32_t) (bases[1] >> 0), 0, (uint32_t) (bases[2] >> 0), 0, (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 32),
				0, (uint32_t) (bases[1] >> 32), 0, (uint32_t) (bases[2] >> 32), 0, (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b[2] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 64),
				0, (uint32_t) (bases[1] >> 64), 0, (uint32_t) (bases[2] >> 64), 0, (uint32_t) (bases[3] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
	b[3] = _mm256_set_epi32(0, (uint32_t) (bases[0] >> 96),
				0, (uint32_t) (bases[1] >> 96), 0, (uint32_t) (bases[2] >> 96), 0, (uint32_t) (bases[3] >> 96));
	one[3] = _mm256_set1_epi64x((uint32_t) (on >> 96));
	p[3] = _mm256_set1_epi64x((uint32_t) (v >> 96));
// p - 1
	avx2_subtract4(m1, p, one);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
	r[3] = b[3];
// MR exponentiation bit per bit
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 0) {
		bit--;
		// Square
		avx2_modsqu4(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx2_modmul4(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[4], mask = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	zero[2] = _mm256_setzero_si256();
	zero[3] = _mm256_setzero_si256();
	if (avx2_cmp_next4(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx2_cmp_next4(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx2_cmp_next4(&mask, r, m1)) {
			return true;
		}
		// square
		avx2_modsqu4(r, p, mmagic);
		// check current result == 1
		if (avx2_neg_cmp_next4(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx2_cmp_next4(&mask, r, m1)) {
		return true;
	}
	return false;
}

#endif				// avx2

// -----------------------------------------------------------------------------------

#define BASES(type, count, mod, var, exp)\
	type var = exp;  \
	while (var >= mod) var -= mod;   \
	montg_bases[count] = var; do { \
	} while (0)

uint32_t montgomery_bases1(uint32_t * montg_bases, uint32_t v, uint64_t count)
{
	// compute x * 2^32 for different values of x

	uint64_t t1 = (uint64_t) 1 << 32;
	t1 %= v;		// 1
	uint64_t t2 = (uint64_t) t1 + t1;	// 2
	while (t2 >= v)
		t2 -= v;
	BASES(uint64_t, 0, v, t3, t1 + t2);
	BASES(uint64_t, 1, v, t5, t3 + t2);
	BASES(uint64_t, 2, v, t7, t5 + t2);
	BASES(uint64_t, 3, v, t11, t5 + t5 + t1);
	if (count <= 4)
		return t1;

	return t1;
}

uint128_t montgomery_bases2(uint64_t * montg_bases, uint64_t v, uint64_t count)
{
	// compute x * 2^64 for different values of x

	uint128_t t1 = (uint128_t) 1 << 64;
	t1 %= v;		// t1

	uint128_t t2 = t1 + t1;	// 2
	while (t2 >= v)
		t2 -= v;

	BASES(uint128_t, 0, v, t3, t1 + t2);
	BASES(uint128_t, 1, v, t5, t3 + t2);
	BASES(uint128_t, 2, v, t7, t5 + t2);
	BASES(uint128_t, 3, v, t11, t5 + t5 + t1);
	if (count <= 4)
		return t1;

	BASES(uint128_t, 4, v, t13, t11 + t2);
	BASES(uint128_t, 5, v, t17, t11 + t5 + t1);
	BASES(uint128_t, 6, v, t19, t17 + t2);
	BASES(uint128_t, 7, v, t23, t11 + t11 + t1);
	if (count <= 8)
		return t1;

	BASES(uint128_t, 8, v, t29, t11 + t17 + t1);
	BASES(uint128_t, 9, v, t31, t29 + t2);
	BASES(uint128_t, 10, v, t37, t17 + t17 + t3);
	BASES(uint128_t, 11, v, t41, t17 + t17 + t7);
	if (count <= 12)
		return t1;

	return t1;
}

uint128_t montgomery_bases3(uint128_t * montg_bases, uint128_t v, uint64_t count)
{
	// compute x * 2^96 for different values of x

	uint128_t t1 = (uint128_t) 1 << 96;
	t1 %= v;		// t1

	uint128_t t2 = t1 + t1;	// 2
	while (t2 >= v)
		t2 -= v;

	BASES(uint128_t, 0, v, t3, t1 + t2);
	BASES(uint128_t, 1, v, t5, t3 + t2);
	BASES(uint128_t, 2, v, t7, t5 + t2);
	BASES(uint128_t, 3, v, t11, t5 + t5 + t1);
	if (count <= 4)
		return t1;

	BASES(uint128_t, 4, v, t13, t11 + t2);
	BASES(uint128_t, 5, v, t17, t11 + t5 + t1);
	BASES(uint128_t, 6, v, t19, t17 + t2);
	BASES(uint128_t, 7, v, t23, t11 + t11 + t1);
	if (count <= 8)
		return t1;

	BASES(uint128_t, 8, v, t29, t11 + t17 + t1);
	BASES(uint128_t, 9, v, t31, t29 + t2);
	BASES(uint128_t, 10, v, t37, t17 + t17 + t3);
	BASES(uint128_t, 11, v, t41, t17 + t17 + t7);
	if (count <= 12)
		return t1;

	BASES(uint128_t, 12, v, t43, t41 + t2);
	BASES(uint128_t, 13, v, t47, t23 + t23 + t1);
	BASES(uint128_t, 14, v, t53, t23 + t23 + t7);
	BASES(uint128_t, 15, v, t59, t23 + t23 + t13);
	if (count <= 16)
		return t1;

	BASES(uint128_t, 16, v, t61, t29 + t31 + t1);
	BASES(uint128_t, 17, v, t67, t29 + t31 + t7);
	BASES(uint128_t, 18, v, t71, t29 + t31 + t11);
	BASES(uint128_t, 19, v, t73, t29 + t31 + t13);
	if (count <= 20)
		return t1;

	return t1;
}

uint128_t montgomery_bases4(uint128_t * montg_bases, uint128_t v, uint64_t count)
{
	// compute x * 2^128 for different values of x

	uint128_t t1 = -v;
	t1 %= v;		// t1

	uint128_t t2 = t1 + t1;	// 2
	while (t2 >= v)
		t2 -= v;

	BASES(uint128_t, 0, v, t3, t1 + t2);
	BASES(uint128_t, 1, v, t5, t3 + t2);
	BASES(uint128_t, 2, v, t7, t5 + t2);
	BASES(uint128_t, 3, v, t11, t5 + t5 + t1);
	if (count <= 4)
		return t1;

	BASES(uint128_t, 4, v, t13, t11 + t2);
	BASES(uint128_t, 5, v, t17, t11 + t5 + t1);
	BASES(uint128_t, 6, v, t19, t17 + t2);
	BASES(uint128_t, 7, v, t23, t11 + t11 + t1);
	if (count <= 8)
		return t1;

	BASES(uint128_t, 8, v, t29, t11 + t17 + t1);
	BASES(uint128_t, 9, v, t31, t29 + t2);
	BASES(uint128_t, 10, v, t37, t17 + t17 + t3);
	BASES(uint128_t, 11, v, t41, t17 + t17 + t7);
	if (count <= 12)
		return t1;

	BASES(uint128_t, 12, v, t43, t41 + t2);
	BASES(uint128_t, 13, v, t47, t23 + t23 + t1);
	BASES(uint128_t, 14, v, t53, t23 + t23 + t7);
	BASES(uint128_t, 15, v, t59, t23 + t23 + t13);
	if (count <= 16)
		return t1;

	BASES(uint128_t, 16, v, t61, t29 + t31 + t1);
	BASES(uint128_t, 17, v, t67, t29 + t31 + t7);
	BASES(uint128_t, 18, v, t71, t29 + t31 + t11);
	BASES(uint128_t, 19, v, t73, t29 + t31 + t13);
	if (count <= 20)
		return t1;

	BASES(uint128_t, 20, v, t79, t29 + t31 + t19);
	BASES(uint128_t, 21, v, t83, t29 + t31 + t23);
	BASES(uint128_t, 22, v, t89, t29 + t31 + t29);
	BASES(uint128_t, 23, v, t97, t29 + t31 + t37);
	if (count <= 24)
		return t1;

	return t1;
}

bool avx2SprpTest(uint64_t v_lo, uint64_t v_hi)
{
	uint128_t v = ((uint128_t) v_hi << 64) + v_lo;

	uint32_t m = montgomeryInverse32((uint32_t) v);

	if (v >> 32 == 0) {
		uint32_t montg_bases[4];
		uint32_t one = montgomery_bases1(montg_bases, v, 4);
		if (!avx2_sprp1(v, m, one, &montg_bases[0]))
			return false;
		return true;
	}
	if (v >> 64 == 0) {
		if (v < 2152302898747ull) {
			uint64_t montg_bases[4];
			uint64_t one = montgomery_bases2(montg_bases, v, 4);
			if (!avx2_sprp2(v, m, one, &montg_bases[0]))
				return false;
			return true;
		} else if (v < 3825123056546413051ull) {
			uint64_t montg_bases[8];
			uint64_t one = montgomery_bases2(montg_bases, v, 8);
			if (!avx2_sprp2(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp2(v, m, one, &montg_bases[4]))
				return false;
			return true;
		} else {
			uint64_t montg_bases[12];
			uint64_t one = montgomery_bases2(montg_bases, v, 12);
			if (!avx2_sprp2(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp2(v, m, one, &montg_bases[4]))
				return false;
			if (!avx2_sprp2(v, m, one, &montg_bases[8]))
				return false;
			return true;
		}
	}
	if (v >> 96 == 0) {
		// if (v < 3317044064679887385961981)
		if (v < ((uint128_t) 0x2BE69ull << 64) + 0x51ADC5B22410A5FDull) {
			uint128_t montg_bases[12];
			uint128_t one = montgomery_bases3(montg_bases, v, 12);
			if (!avx2_sprp3(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp3(v, m, one, &montg_bases[4]))
				return false;
			if (!avx2_sprp3(v, m, one, &montg_bases[8]))
				return false;
			return true;
		} else {
			uint128_t montg_bases[16];
			uint128_t one = montgomery_bases3(montg_bases, v, 16);
			if (!avx2_sprp3(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp3(v, m, one, &montg_bases[4]))
				return false;
			if (!avx2_sprp3(v, m, one, &montg_bases[8]))
				return false;
			if (!avx2_sprp3(v, m, one, &montg_bases[12]))
				return false;
			return true;
		}
	} else {
		if (v >> 112 == 0) {
			uint128_t montg_bases[16];
			uint128_t one = montgomery_bases4(montg_bases, v, 16);
			if (!avx2_sprp4(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[4]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[8]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[12]))
				return false;
			return true;
		} else {
			uint128_t montg_bases[20];
			uint128_t one = montgomery_bases4(montg_bases, v, 20);
			if (!avx2_sprp4(v, m, one, &montg_bases[0]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[4]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[8]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[12]))
				return false;
			if (!avx2_sprp4(v, m, one, &montg_bases[16]))
				return false;
			return true;
		}
	}

	return false;

	/*
	   bases[0] = _mm256_set_epi64(4130806001517ull, 149795463772692060ull, 186635894390467037ull, 3967304179347715805ull);
	   bases[1] = _mm256_set_epi64(4130806001517ull >> 32, 149795463772692060ull >> 32, 186635894390467037ull >>32, 3967304179347715805ull >> 32);
	 */
}
