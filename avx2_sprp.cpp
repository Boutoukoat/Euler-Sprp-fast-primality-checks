
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <x86intrin.h>

#include "avx2_sprp.h"
#include "m_reg.h"
#include "m128_utils.h"
#include "m64_utils.h"

#if defined(__AVX2__)

// -----------------------------------------------------------------------------------
//
//Generated code starts here
//
// -----------------------------------------------------------------------------------

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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx2_neg_cmp_next1(__m256i * mask, __m256i ma[1], __m256i mb[1])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
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

static
    inline __attribute__((always_inline))
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
	__m256i p[1], r[1], one[1], m1[1], b[1], bbb[1];
	uint64_t bit, k;
	uint32_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				 (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
// p - 1
	avx2_subtract1(m1, p, one);
// windowing bbb = b^3 mod p
	bbb[0] = b[0];
	avx2_modsqu1(bbb, p, mmagic);
	avx2_modmul1(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
// exponentiation bit per bit (with sliding window size 2)
	k = my_ctz32(v - 1);
	s = v >> k;
	bit = 31 - my_clz32(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx2_modsqu1(r, p, mmagic);
			avx2_modsqu1(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx2_modsqu1(r, p, mmagic);
			break;
		case 2:
			avx2_modsqu1(r, p, mmagic);
			avx2_modmul1(r, b, p, mmagic);
			avx2_modsqu1(r, p, mmagic);
			break;
		case 3:
			avx2_modsqu1(r, p, mmagic);
			avx2_modsqu1(r, p, mmagic);
			avx2_modmul1(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
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
// final iteration square per square
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
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

static
    inline __attribute__((always_inline))
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
	__m256i p[2], r[2], one[2], m1[2], b[2], bbb[2];
	uint64_t bit, k;
	uint64_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				 (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				 (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
// p - 1
	avx2_subtract2(m1, p, one);
// windowing bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	avx2_modsqu2(bbb, p, mmagic);
	avx2_modmul2(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
// exponentiation bit per bit (with sliding window size 2)
	k = my_ctz64(v - 1);
	s = v >> k;
	bit = 63 - my_clz64(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx2_modsqu2(r, p, mmagic);
			avx2_modsqu2(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx2_modsqu2(r, p, mmagic);
			break;
		case 2:
			avx2_modsqu2(r, p, mmagic);
			avx2_modmul2(r, b, p, mmagic);
			avx2_modsqu2(r, p, mmagic);
			break;
		case 3:
			avx2_modsqu2(r, p, mmagic);
			avx2_modsqu2(r, p, mmagic);
			avx2_modmul2(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
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
// final iteration square per square
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
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

static
    inline __attribute__((always_inline))
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
	__m256i p[3], r[3], one[3], m1[3], b[3], bbb[3];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				 (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				 (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b[2] = _mm256_set_epi64x((uint32_t) (bases[0] >> 64),
				 (uint32_t) (bases[1] >> 64), (uint32_t) (bases[2] >> 64), (uint32_t) (bases[3] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
// p - 1
	avx2_subtract3(m1, p, one);
// windowing bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	bbb[2] = b[2];
	avx2_modsqu3(bbb, p, mmagic);
	avx2_modmul3(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
// exponentiation bit per bit (with sliding window size 2)
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx2_modsqu3(r, p, mmagic);
			avx2_modsqu3(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx2_modsqu3(r, p, mmagic);
			break;
		case 2:
			avx2_modsqu3(r, p, mmagic);
			avx2_modmul3(r, b, p, mmagic);
			avx2_modsqu3(r, p, mmagic);
			break;
		case 3:
			avx2_modsqu3(r, p, mmagic);
			avx2_modsqu3(r, p, mmagic);
			avx2_modmul3(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
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
// final iteration square per square
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[3], mb[3]);
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
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[3], mb[3]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
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

static
    inline __attribute__((always_inline))
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
	__m256i p[4], r[4], one[4], m1[4], b[4], bbb[4];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				 (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				 (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b[2] = _mm256_set_epi64x((uint32_t) (bases[0] >> 64),
				 (uint32_t) (bases[1] >> 64), (uint32_t) (bases[2] >> 64), (uint32_t) (bases[3] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
	b[3] = _mm256_set_epi64x((uint32_t) (bases[0] >> 96),
				 (uint32_t) (bases[1] >> 96), (uint32_t) (bases[2] >> 96), (uint32_t) (bases[3] >> 96));
	one[3] = _mm256_set1_epi64x((uint32_t) (on >> 96));
	p[3] = _mm256_set1_epi64x((uint32_t) (v >> 96));
// p - 1
	avx2_subtract4(m1, p, one);
// windowing bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	bbb[2] = b[2];
	bbb[3] = b[3];
	avx2_modsqu4(bbb, p, mmagic);
	avx2_modmul4(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
	r[3] = b[3];
// exponentiation bit per bit (with sliding window size 2)
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx2_modsqu4(r, p, mmagic);
			avx2_modsqu4(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx2_modsqu4(r, p, mmagic);
			break;
		case 2:
			avx2_modsqu4(r, p, mmagic);
			avx2_modmul4(r, b, p, mmagic);
			avx2_modsqu4(r, p, mmagic);
			break;
		case 3:
			avx2_modsqu4(r, p, mmagic);
			avx2_modsqu4(r, p, mmagic);
			avx2_modmul4(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
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
// final iteration square per square
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
//
//Generated code ends here
//
// -----------------------------------------------------------------------------------


#if defined(__AVX2__)

// -----------------------------------------------------------------------------------
//
//Generated code starts here
//
// -----------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------

// if (mr1 >= mp) mr1 -= mp
// if (mr2 >= mp) mr2 -= mp
static inline __attribute__((always_inline))
void avx22_modsubtract1(__m256i * mr1, __m256i * mr2, __m256i mp[1])
{
	__m256i t1[1], t2[1], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx22_modsubtract21(__m256i * mr1, __m256i mc1, __m256i * mr2, __m256i mc2, __m256i mp[1])
{
	__m256i t1[1], t2[1], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mc1);
	cs2 = _mm256_add_epi64(cs2, mc2);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
}

// -----------------------------------------------------------------------------------

// mr = ma - mb
static inline __attribute__((always_inline))
void avx22_subtract1(__m256i * mr, __m256i ma[1], __m256i mb[1])
{
	mr[0] = _mm256_sub_epi64(ma[0], mb[0]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx22_cmp_next1(__m256i * mask, __m256i ma[1], __m256i mb[1])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx22_neg_cmp_next1(__m256i * mask, __m256i ma[1], __m256i mb[1])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modmul1(__m256i * mr1, __m256i mb1[1], __m256i * mr2, __m256i mb2[1], __m256i mp[1], __m256i mmagic)
{
	__m256i t1[3];
	__m256i t2[3];
	__m256i cm1, c1, cs1, d1, ds1, a1;
	__m256i cm2, c2, cs2, d2, ds2, a2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #1
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #3
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #4
	t1[1] = _mm256_srli_epi64(cs1, 32);
	t2[1] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #5
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #6
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, t1[1]);
	ds2 = _mm256_add_epi64(d2, t2[1]);
	mr1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	ds1 = _mm256_srli_epi64(ds1, 32);
	ds2 = _mm256_srli_epi64(ds2, 32);
	avx22_modsubtract21(mr1, ds1, mr2, ds2, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modsqu1(__m256i * mr1, __m256i * mr2, __m256i mp[1], __m256i mmagic)
{
	__m256i t1[3];
	__m256i t2[3];
	__m256i a1, cm1, c1, ca1, cb1, cs1, ct1, cd1;
	__m256i a2, cm2, c2, ca2, cb2, cs2, ct2, cd2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, a1);	// #1
	cs2 = _mm256_mul_epu32(a2, a2);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #3
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #4
	t1[1] = _mm256_srli_epi64(cs1, 32);
	t2[1] = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #5
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #6
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[1]);
	cs2 = _mm256_add_epi64(c2, t2[1]);
	mr1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_srli_epi64(cs1, 32);
	cs2 = _mm256_srli_epi64(cs2, 32);
	avx22_modsubtract21(mr1, cs1, mr2, cs2, mp);
}

// -----------------------------------------------------------------------------------

bool avx22_sprp1(uint32_t v, uint32_t mm, uint32_t on, uint32_t * bases)
{
	__m256i p[1], r1[1], r2[1], one[1], m1[1], b1[1], b2[1], bbb1[1], bbb2[1];
	uint64_t bit, k;
	uint32_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b1[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				  (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	b2[0] = _mm256_set_epi64x((uint32_t) (bases[4] >> 0),
				  (uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
// p - 1
	avx22_subtract1(m1, p, one);
// first value
	r1[0] = b1[0];
	r2[0] = b2[0];
// Windowing
	bbb1[0] = b1[0];
	bbb2[0] = b2[0];
	avx22_modsqu1(bbb1, bbb2, p, mmagic);
	avx22_modmul1(bbb1, b1, bbb2, b2, p, mmagic);
// MR exponentiation bit per bit
	k = my_ctz32(v - 1);
	s = v >> k;
	bit = 31 - my_clz32(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx22_modsqu1(r1, r2, p, mmagic);
			avx22_modsqu1(r1, r2, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx22_modsqu1(r1, r2, p, mmagic);
			break;
		case 2:
			avx22_modsqu1(r1, r2, p, mmagic);
			avx22_modmul1(r1, b1, r2, b2, p, mmagic);
			avx22_modsqu1(r1, r2, p, mmagic);
			break;
		case 3:
			avx22_modsqu1(r1, r2, p, mmagic);
			avx22_modsqu1(r1, r2, p, mmagic);
			avx22_modmul1(r1, bbb1, r2, bbb2, p, mmagic);
			break;
		}
	}
	while (bit > 0) {
		bit--;
		// Square
		avx22_modsqu1(r1, r2, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx22_modmul1(r1, b1, r2, b2, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[1];
	__m256i mask1 = _mm256_setzero_si256();
	__m256i mask2 = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	if (avx22_cmp_next1(&mask1, b1, zero)) {
		return true;
	}
	if (avx22_cmp_next1(&mask2, b2, zero)) {
		return true;
	}
// check current result == 1
	if (avx22_cmp_next1(&mask1, r1, one)) {
		return true;
	}
	if (avx22_cmp_next1(&mask2, r2, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx22_cmp_next1(&mask1, r1, m1)) {
			return true;
		}
		if (avx22_cmp_next1(&mask2, r2, m1)) {
			return true;
		}
		// square
		avx22_modsqu1(r1, r2, p, mmagic);
		// check current result == 1
		if (avx22_neg_cmp_next1(&mask1, r1, one)) {
			// non-trivial quadratic residue found
			return false;
		}
		if (avx22_neg_cmp_next1(&mask2, r2, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx22_cmp_next1(&mask1, r1, m1)) {
		return true;
	}
	if (avx22_cmp_next1(&mask2, r2, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

// if (mr1 >= mp) mr1 -= mp
// if (mr2 >= mp) mr2 -= mp
static inline __attribute__((always_inline))
void avx22_modsubtract2(__m256i * mr1, __m256i * mr2, __m256i mp[2])
{
	__m256i t1[2], t2[2], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx22_modsubtract32(__m256i * mr1, __m256i mc1, __m256i * mr2, __m256i mc2, __m256i mp[2])
{
	__m256i t1[2], t2[2], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mc1);
	cs2 = _mm256_add_epi64(cs2, mc2);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
}

// -----------------------------------------------------------------------------------

// mr = ma - mb
static inline __attribute__((always_inline))
void avx22_subtract2(__m256i * mr, __m256i ma[2], __m256i mb[2])
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
bool avx22_cmp_next2(__m256i * mask, __m256i ma[2], __m256i mb[2])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx22_neg_cmp_next2(__m256i * mask, __m256i ma[2], __m256i mb[2])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modmul2(__m256i * mr1, __m256i mb1[2], __m256i * mr2, __m256i mb2[2], __m256i mp[2], __m256i mmagic)
{
	__m256i t1[4];
	__m256i t2[4];
	__m256i cm1, c1, cs1, d1, ds1, a1;
	__m256i cm2, c2, cs2, d2, ds2, a2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #1
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #3
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #4
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #5
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #6
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #7
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #8
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#9
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#10
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, c1);
	ds2 = _mm256_add_epi64(d2, c2);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	t1[2] = _mm256_srli_epi64(ds1, 32);
	t2[2] = _mm256_srli_epi64(ds2, 32);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #11
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #12
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #13
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #14
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #15
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #16
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #17
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #18
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#19
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#20
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	mr1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[2]);
	cs2 = _mm256_add_epi64(c2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[3] = _mm256_srli_epi64(cs1, 32);
	t2[3] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[2]);
	ds2 = _mm256_add_epi64(d2, t2[2]);
	mr1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, t1[3]);
	ds2 = _mm256_add_epi64(d2, t2[3]);

	avx22_modsubtract32(mr1, ds1, mr2, ds2, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modsqu2(__m256i * mr1, __m256i * mr2, __m256i mp[2], __m256i mmagic)
{
	__m256i t1[4];
	__m256i t2[4];
	__m256i a1, cm1, c1, ca1, cb1, cs1, ct1, cd1;
	__m256i a2, cm2, c2, ca2, cb2, cs2, ct2, cd2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, a1);	// #1
	cs2 = _mm256_mul_epu32(a2, a2);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[1]);	// #3
	ct2 = _mm256_mul_epu32(a2, mr2[1]);	// #4
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[3] = _mm256_srli_epi64(cs1, 32);
	t2[3] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #5
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #6
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #7
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #8
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #9
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #10
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[2]);
	cs2 = _mm256_add_epi64(c2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[2] = _mm256_add_epi64(c1, t1[3]);
	t2[2] = _mm256_add_epi64(c2, t2[3]);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, a1);	// #11
	cs2 = _mm256_mul_epu32(a2, a2);	// #12
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ca1, t1[2]);
	cs2 = _mm256_add_epi64(ca2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[3] = _mm256_srli_epi64(cs1, 32);
	t2[3] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #13
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #14
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #15
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #16
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #17
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #18
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	mr1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[2]);
	cs2 = _mm256_add_epi64(c2, t2[2]);
	mr1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);

	avx22_modsubtract32(mr1, cs1, mr2, cs2, mp);
}

// -----------------------------------------------------------------------------------

bool avx22_sprp2(uint64_t v, uint32_t mm, uint64_t on, uint64_t * bases)
{
	__m256i p[2], r1[2], r2[2], one[2], m1[2], b1[2], b2[2], bbb1[2], bbb2[2];
	uint64_t bit, k;
	uint64_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b1[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				  (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	b2[0] = _mm256_set_epi64x((uint32_t) (bases[4] >> 0),
				  (uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b1[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				  (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	b2[1] = _mm256_set_epi64x((uint32_t) (bases[4] >> 32),
				  (uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
// p - 1
	avx22_subtract2(m1, p, one);
// first value
	r1[0] = b1[0];
	r2[0] = b2[0];
	r1[1] = b1[1];
	r2[1] = b2[1];
// Windowing
	bbb1[0] = b1[0];
	bbb2[0] = b2[0];
	bbb1[1] = b1[1];
	bbb2[1] = b2[1];
	avx22_modsqu2(bbb1, bbb2, p, mmagic);
	avx22_modmul2(bbb1, b1, bbb2, b2, p, mmagic);
// MR exponentiation bit per bit
	k = my_ctz64(v - 1);
	s = v >> k;
	bit = 63 - my_clz64(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx22_modsqu2(r1, r2, p, mmagic);
			avx22_modsqu2(r1, r2, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx22_modsqu2(r1, r2, p, mmagic);
			break;
		case 2:
			avx22_modsqu2(r1, r2, p, mmagic);
			avx22_modmul2(r1, b1, r2, b2, p, mmagic);
			avx22_modsqu2(r1, r2, p, mmagic);
			break;
		case 3:
			avx22_modsqu2(r1, r2, p, mmagic);
			avx22_modsqu2(r1, r2, p, mmagic);
			avx22_modmul2(r1, bbb1, r2, bbb2, p, mmagic);
			break;
		}
	}
	while (bit > 0) {
		bit--;
		// Square
		avx22_modsqu2(r1, r2, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx22_modmul2(r1, b1, r2, b2, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[2];
	__m256i mask1 = _mm256_setzero_si256();
	__m256i mask2 = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	if (avx22_cmp_next2(&mask1, b1, zero)) {
		return true;
	}
	if (avx22_cmp_next2(&mask2, b2, zero)) {
		return true;
	}
// check current result == 1
	if (avx22_cmp_next2(&mask1, r1, one)) {
		return true;
	}
	if (avx22_cmp_next2(&mask2, r2, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx22_cmp_next2(&mask1, r1, m1)) {
			return true;
		}
		if (avx22_cmp_next2(&mask2, r2, m1)) {
			return true;
		}
		// square
		avx22_modsqu2(r1, r2, p, mmagic);
		// check current result == 1
		if (avx22_neg_cmp_next2(&mask1, r1, one)) {
			// non-trivial quadratic residue found
			return false;
		}
		if (avx22_neg_cmp_next2(&mask2, r2, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx22_cmp_next2(&mask1, r1, m1)) {
		return true;
	}
	if (avx22_cmp_next2(&mask2, r2, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

// if (mr1 >= mp) mr1 -= mp
// if (mr2 >= mp) mr2 -= mp
static inline __attribute__((always_inline))
void avx22_modsubtract3(__m256i * mr1, __m256i * mr2, __m256i mp[3])
{
	__m256i t1[3], t2[3], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[2]);
	cs2 = _mm256_add_epi64(cs2, mr2[2]);
	cs1 = _mm256_sub_epi64(cs1, mp[2]);
	cs2 = _mm256_sub_epi64(cs2, mp[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
	mr1[2] = _mm256_blendv_epi8(t1[2], mr1[2], cs1);
	mr2[2] = _mm256_blendv_epi8(t2[2], mr2[2], cs2);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx22_modsubtract43(__m256i * mr1, __m256i mc1, __m256i * mr2, __m256i mc2, __m256i mp[3])
{
	__m256i t1[3], t2[3], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[2]);
	cs2 = _mm256_add_epi64(cs2, mr2[2]);
	cs1 = _mm256_sub_epi64(cs1, mp[2]);
	cs2 = _mm256_sub_epi64(cs2, mp[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mc1);
	cs2 = _mm256_add_epi64(cs2, mc2);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
	mr1[2] = _mm256_blendv_epi8(t1[2], mr1[2], cs1);
	mr2[2] = _mm256_blendv_epi8(t2[2], mr2[2], cs2);
}

// -----------------------------------------------------------------------------------

// mr = ma - mb
static inline __attribute__((always_inline))
void avx22_subtract3(__m256i * mr, __m256i ma[3], __m256i mb[3])
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
bool avx22_cmp_next3(__m256i * mask, __m256i ma[3], __m256i mb[3])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx22_neg_cmp_next3(__m256i * mask, __m256i ma[3], __m256i mb[3])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modmul3(__m256i * mr1, __m256i mb1[3], __m256i * mr2, __m256i mb2[3], __m256i mp[3], __m256i mmagic)
{
	__m256i t1[5];
	__m256i t2[5];
	__m256i cm1, c1, cs1, d1, ds1, a1;
	__m256i cm2, c2, cs2, d2, ds2, a2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #1
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #3
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #4
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #5
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #6
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #7
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #8
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#9
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#10
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #11
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #12
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#13
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#14
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, c1);
	ds2 = _mm256_add_epi64(d2, c2);
	t1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	t1[3] = _mm256_srli_epi64(ds1, 32);
	t2[3] = _mm256_srli_epi64(ds2, 32);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #15
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #16
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #17
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #18
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #19
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #20
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #21
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #22
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#23
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#24
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #25
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #26
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#27
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#28
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(cs1, 32);
	t2[4] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[3]);
	ds2 = _mm256_add_epi64(d2, t2[3]);
	t1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	t1[3] = _mm256_add_epi64(d1, t1[4]);
	t2[3] = _mm256_add_epi64(d2, t2[4]);

	a1 = mr1[2];
	a2 = mr2[2];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #29
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #30
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #31
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #32
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #33
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #34
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #35
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #36
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#37
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#38
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	mr1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #39
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #40
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#41
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#42
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	mr1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(cs1, 32);
	t2[4] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[3]);
	ds2 = _mm256_add_epi64(d2, t2[3]);
	mr1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, t1[4]);
	ds2 = _mm256_add_epi64(d2, t2[4]);

	avx22_modsubtract43(mr1, ds1, mr2, ds2, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modsqu3(__m256i * mr1, __m256i * mr2, __m256i mp[3], __m256i mmagic)
{
	__m256i t1[5];
	__m256i t2[5];
	__m256i a1, cm1, c1, ca1, cb1, cs1, ct1, cd1;
	__m256i a2, cm2, c2, ca2, cb2, cs2, ct2, cd2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, a1);	// #1
	cs2 = _mm256_mul_epu32(a2, a2);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[1]);	// #3
	ct2 = _mm256_mul_epu32(a2, mr2[1]);	// #4
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[2]);	// #5
	ct2 = _mm256_mul_epu32(a2, mr2[2]);	// #6
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	cs1 = _mm256_add_epi64(cs1, cb1);
	cs2 = _mm256_add_epi64(cs2, cb2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(cs1, 32);
	t2[4] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #7
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #8
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #9
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #10
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #11
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #12
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #13
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #14
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[3] = _mm256_add_epi64(c1, t1[4]);
	t2[3] = _mm256_add_epi64(c2, t2[4]);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, a1);	// #15
	cs2 = _mm256_mul_epu32(a2, a2);	// #16
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[2]);	// #17
	ct2 = _mm256_mul_epu32(a2, mr2[2]);	// #18
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(cs1, 32);
	t2[4] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #19
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #20
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #21
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #22
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #23
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #24
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #25
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #26
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[3] = _mm256_add_epi64(c1, t1[4]);
	t2[3] = _mm256_add_epi64(c2, t2[4]);

	a1 = mr1[2];
	a2 = mr2[2];
	cs1 = _mm256_mul_epu32(a1, a1);	// #27
	cs2 = _mm256_mul_epu32(a2, a2);	// #28
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ca1, t1[3]);
	cs2 = _mm256_add_epi64(ca2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(cs1, 32);
	t2[4] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #29
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #30
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #31
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #32
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #33
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #34
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	mr1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #35
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #36
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	mr1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[3]);
	cs2 = _mm256_add_epi64(c2, t2[3]);
	mr1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);

	avx22_modsubtract43(mr1, cs1, mr2, cs2, mp);
}

// -----------------------------------------------------------------------------------

bool avx22_sprp3(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m256i p[3], r1[3], r2[3], one[3], m1[3], b1[3], b2[3], bbb1[3], bbb2[3];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b1[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				  (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	b2[0] = _mm256_set_epi64x((uint32_t) (bases[4] >> 0),
				  (uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b1[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				  (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	b2[1] = _mm256_set_epi64x((uint32_t) (bases[4] >> 32),
				  (uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b1[2] = _mm256_set_epi64x((uint32_t) (bases[0] >> 64),
				  (uint32_t) (bases[1] >> 64), (uint32_t) (bases[2] >> 64), (uint32_t) (bases[3] >> 64));
	b2[2] = _mm256_set_epi64x((uint32_t) (bases[4] >> 64),
				  (uint32_t) (bases[5] >> 64), (uint32_t) (bases[6] >> 64), (uint32_t) (bases[7] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
// p - 1
	avx22_subtract3(m1, p, one);
// first value
	r1[0] = b1[0];
	r2[0] = b2[0];
	r1[1] = b1[1];
	r2[1] = b2[1];
	r1[2] = b1[2];
	r2[2] = b2[2];
// Windowing
	bbb1[0] = b1[0];
	bbb2[0] = b2[0];
	bbb1[1] = b1[1];
	bbb2[1] = b2[1];
	bbb1[2] = b1[2];
	bbb2[2] = b2[2];
	avx22_modsqu3(bbb1, bbb2, p, mmagic);
	avx22_modmul3(bbb1, b1, bbb2, b2, p, mmagic);
// MR exponentiation bit per bit
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx22_modsqu3(r1, r2, p, mmagic);
			avx22_modsqu3(r1, r2, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx22_modsqu3(r1, r2, p, mmagic);
			break;
		case 2:
			avx22_modsqu3(r1, r2, p, mmagic);
			avx22_modmul3(r1, b1, r2, b2, p, mmagic);
			avx22_modsqu3(r1, r2, p, mmagic);
			break;
		case 3:
			avx22_modsqu3(r1, r2, p, mmagic);
			avx22_modsqu3(r1, r2, p, mmagic);
			avx22_modmul3(r1, bbb1, r2, bbb2, p, mmagic);
			break;
		}
	}
	while (bit > 0) {
		bit--;
		// Square
		avx22_modsqu3(r1, r2, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx22_modmul3(r1, b1, r2, b2, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[3];
	__m256i mask1 = _mm256_setzero_si256();
	__m256i mask2 = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	zero[2] = _mm256_setzero_si256();
	if (avx22_cmp_next3(&mask1, b1, zero)) {
		return true;
	}
	if (avx22_cmp_next3(&mask2, b2, zero)) {
		return true;
	}
// check current result == 1
	if (avx22_cmp_next3(&mask1, r1, one)) {
		return true;
	}
	if (avx22_cmp_next3(&mask2, r2, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx22_cmp_next3(&mask1, r1, m1)) {
			return true;
		}
		if (avx22_cmp_next3(&mask2, r2, m1)) {
			return true;
		}
		// square
		avx22_modsqu3(r1, r2, p, mmagic);
		// check current result == 1
		if (avx22_neg_cmp_next3(&mask1, r1, one)) {
			// non-trivial quadratic residue found
			return false;
		}
		if (avx22_neg_cmp_next3(&mask2, r2, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx22_cmp_next3(&mask1, r1, m1)) {
		return true;
	}
	if (avx22_cmp_next3(&mask2, r2, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

// if (mr1 >= mp) mr1 -= mp
// if (mr2 >= mp) mr2 -= mp
static inline __attribute__((always_inline))
void avx22_modsubtract4(__m256i * mr1, __m256i * mr2, __m256i mp[4])
{
	__m256i t1[4], t2[4], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[2]);
	cs2 = _mm256_add_epi64(cs2, mr2[2]);
	cs1 = _mm256_sub_epi64(cs1, mp[2]);
	cs2 = _mm256_sub_epi64(cs2, mp[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[3]);
	cs2 = _mm256_add_epi64(cs2, mr2[3]);
	cs1 = _mm256_sub_epi64(cs1, mp[3]);
	cs2 = _mm256_sub_epi64(cs2, mp[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
	mr1[2] = _mm256_blendv_epi8(t1[2], mr1[2], cs1);
	mr2[2] = _mm256_blendv_epi8(t2[2], mr2[2], cs2);
	mr1[3] = _mm256_blendv_epi8(t1[3], mr1[3], cs1);
	mr2[3] = _mm256_blendv_epi8(t2[3], mr2[3], cs2);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx22_modsubtract54(__m256i * mr1, __m256i mc1, __m256i * mr2, __m256i mc2, __m256i mp[4])
{
	__m256i t1[4], t2[4], cs1, cs2, zero = _mm256_setzero_si256();
	cs1 = _mm256_sub_epi64(mr1[0], mp[0]);
	cs2 = _mm256_sub_epi64(mr2[0], mp[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[1]);
	cs2 = _mm256_add_epi64(cs2, mr2[1]);
	cs1 = _mm256_sub_epi64(cs1, mp[1]);
	cs2 = _mm256_sub_epi64(cs2, mp[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[2]);
	cs2 = _mm256_add_epi64(cs2, mr2[2]);
	cs1 = _mm256_sub_epi64(cs1, mp[2]);
	cs2 = _mm256_sub_epi64(cs2, mp[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mr1[3]);
	cs2 = _mm256_add_epi64(cs2, mr2[3]);
	cs1 = _mm256_sub_epi64(cs1, mp[3]);
	cs2 = _mm256_sub_epi64(cs2, mp[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	cs1 = _mm256_add_epi64(cs1, mc1);
	cs2 = _mm256_add_epi64(cs2, mc2);
	cs1 = _mm256_shuffle_epi32(cs1, _MM_SHUFFLE(3, 3, 1, 1));
	cs2 = _mm256_shuffle_epi32(cs2, _MM_SHUFFLE(3, 3, 1, 1));
	mr1[0] = _mm256_blendv_epi8(t1[0], mr1[0], cs1);
	mr2[0] = _mm256_blendv_epi8(t2[0], mr2[0], cs2);
	mr1[1] = _mm256_blendv_epi8(t1[1], mr1[1], cs1);
	mr2[1] = _mm256_blendv_epi8(t2[1], mr2[1], cs2);
	mr1[2] = _mm256_blendv_epi8(t1[2], mr1[2], cs1);
	mr2[2] = _mm256_blendv_epi8(t2[2], mr2[2], cs2);
	mr1[3] = _mm256_blendv_epi8(t1[3], mr1[3], cs1);
	mr2[3] = _mm256_blendv_epi8(t2[3], mr2[3], cs2);
}

// -----------------------------------------------------------------------------------

// mr = ma - mb
static inline __attribute__((always_inline))
void avx22_subtract4(__m256i * mr, __m256i ma[4], __m256i mb[4])
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
bool avx22_cmp_next4(__m256i * mask, __m256i ma[4], __m256i mb[4])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[3], mb[3]);
	t0 = _mm256_and_si256(t0, t1);
	*mask = _mm256_or_si256(t0, *mask);
	unsigned f = _mm256_movemask_epi8(*mask);
	return f == (unsigned)0xffffffff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx22_neg_cmp_next4(__m256i * mask, __m256i ma[4], __m256i mb[4])
{
	__m256i t0, t1;
	t0 = _mm256_cmpeq_epi64(ma[0], mb[0]);
	t1 = _mm256_cmpeq_epi64(ma[1], mb[1]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[2], mb[2]);
	t0 = _mm256_and_si256(t0, t1);
	t1 = _mm256_cmpeq_epi64(ma[3], mb[3]);
	t0 = _mm256_and_si256(t0, t1);
	unsigned f = _mm256_testc_si256(*mask, t0);
	return f == 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modmul4(__m256i * mr1, __m256i mb1[4], __m256i * mr2, __m256i mb2[4], __m256i mp[4], __m256i mmagic)
{
	__m256i t1[6];
	__m256i t2[6];
	__m256i cm1, c1, cs1, d1, ds1, a1;
	__m256i cm2, c2, cs2, d2, ds2, a2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #1
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #3
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #4
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #5
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #6
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #7
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #8
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#9
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#10
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #11
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #12
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#13
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#14
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[3]);	// #15
	cs2 = _mm256_mul_epu32(a2, mb2[3]);	// #16
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[3]);	//#17
	ds2 = _mm256_mul_epu32(cm2, mp[3]);	//#18
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[3]);
	ds2 = _mm256_add_epi64(ds2, t2[3]);
	t1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, c1);
	ds2 = _mm256_add_epi64(d2, c2);
	t1[3] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(ds2, zero, 0xaa);
	t1[4] = _mm256_srli_epi64(ds1, 32);
	t2[4] = _mm256_srli_epi64(ds2, 32);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #19
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #20
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #21
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #22
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #23
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #24
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #25
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #26
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#27
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#28
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #29
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #30
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#31
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#32
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[3]);	// #33
	cs2 = _mm256_mul_epu32(a2, mb2[3]);	// #34
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[3]);	//#35
	ds2 = _mm256_mul_epu32(cm2, mp[3]);	//#36
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[3]);
	ds2 = _mm256_add_epi64(ds2, t2[3]);
	t1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[4]);
	ds2 = _mm256_add_epi64(d2, t2[4]);
	t1[3] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	t1[4] = _mm256_add_epi64(d1, t1[5]);
	t2[4] = _mm256_add_epi64(d2, t2[5]);

	a1 = mr1[2];
	a2 = mr2[2];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #37
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #38
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #39
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #40
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #41
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #42
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #43
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #44
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#45
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#46
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	t1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #47
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #48
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#49
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#50
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	t1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[3]);	// #51
	cs2 = _mm256_mul_epu32(a2, mb2[3]);	// #52
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[3]);	//#53
	ds2 = _mm256_mul_epu32(cm2, mp[3]);	//#54
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[3]);
	ds2 = _mm256_add_epi64(ds2, t2[3]);
	t1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[4]);
	ds2 = _mm256_add_epi64(d2, t2[4]);
	t1[3] = _mm256_blend_epi32(ds1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	t1[4] = _mm256_add_epi64(d1, t1[5]);
	t2[4] = _mm256_add_epi64(d2, t2[5]);

	a1 = mr1[3];
	a2 = mr2[3];
	cs1 = _mm256_mul_epu32(a1, mb1[0]);	// #55
	cs2 = _mm256_mul_epu32(a2, mb2[0]);	// #56
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #57
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #58
	ds1 = _mm256_mul_epu32(cm1, mp[0]);	// #59
	ds2 = _mm256_mul_epu32(cm2, mp[0]);	// #60
	ds1 = _mm256_add_epi64(ds1, t1[0]);
	ds2 = _mm256_add_epi64(ds2, t2[0]);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[1]);	// #61
	cs2 = _mm256_mul_epu32(a2, mb2[1]);	// #62
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[1]);	//#63
	ds2 = _mm256_mul_epu32(cm2, mp[1]);	//#64
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[1]);
	ds2 = _mm256_add_epi64(ds2, t2[1]);
	mr1[0] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[2]);	// #65
	cs2 = _mm256_mul_epu32(a2, mb2[2]);	// #66
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[2]);	//#67
	ds2 = _mm256_mul_epu32(cm2, mp[2]);	//#68
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[2]);
	ds2 = _mm256_add_epi64(ds2, t2[2]);
	mr1[1] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_mul_epu32(a1, mb1[3]);	// #69
	cs2 = _mm256_mul_epu32(a2, mb2[3]);	// #70
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_mul_epu32(cm1, mp[3]);	//#71
	ds2 = _mm256_mul_epu32(cm2, mp[3]);	//#72
	ds1 = _mm256_add_epi64(ds1, d1);
	ds2 = _mm256_add_epi64(ds2, d2);
	ds1 = _mm256_add_epi64(ds1, t1[3]);
	ds2 = _mm256_add_epi64(ds2, t2[3]);
	mr1[2] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[2] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);
	ds1 = _mm256_add_epi64(d1, t1[4]);
	ds2 = _mm256_add_epi64(d2, t2[4]);
	mr1[3] = _mm256_blend_epi32(ds1, zero, 0xaa);
	mr2[3] = _mm256_blend_epi32(ds2, zero, 0xaa);
	d1 = _mm256_srli_epi64(ds1, 32);
	d2 = _mm256_srli_epi64(ds2, 32);
	ds1 = _mm256_add_epi64(d1, t1[5]);
	ds2 = _mm256_add_epi64(d2, t2[5]);

	avx22_modsubtract54(mr1, ds1, mr2, ds2, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx22_modsqu4(__m256i * mr1, __m256i * mr2, __m256i mp[4], __m256i mmagic)
{
	__m256i t1[6];
	__m256i t2[6];
	__m256i a1, cm1, c1, ca1, cb1, cs1, ct1, cd1;
	__m256i a2, cm2, c2, ca2, cb2, cs2, ct2, cd2;
	__m256i zero = _mm256_setzero_si256();

	a1 = mr1[0];
	a2 = mr2[0];
	cs1 = _mm256_mul_epu32(a1, a1);	// #1
	cs2 = _mm256_mul_epu32(a2, a2);	// #2
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[1]);	// #3
	ct2 = _mm256_mul_epu32(a2, mr2[1]);	// #4
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[2]);	// #5
	ct2 = _mm256_mul_epu32(a2, mr2[2]);	// #6
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	cs1 = _mm256_add_epi64(cs1, cb1);
	cs2 = _mm256_add_epi64(cs2, cb2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[3]);	// #7
	ct2 = _mm256_mul_epu32(a2, mr2[3]);	// #8
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	cs1 = _mm256_add_epi64(cs1, cb1);
	cs2 = _mm256_add_epi64(cs2, cb2);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #9
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #10
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #11
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #12
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #13
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #14
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #15
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #16
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[3]);	// #17
	cs2 = _mm256_mul_epu32(cm2, mp[3]);	// #18
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[4] = _mm256_add_epi64(c1, t1[5]);
	t2[4] = _mm256_add_epi64(c2, t2[5]);

	a1 = mr1[1];
	a2 = mr2[1];
	cs1 = _mm256_mul_epu32(a1, a1);	// #19
	cs2 = _mm256_mul_epu32(a2, a2);	// #20
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[2]);	// #21
	ct2 = _mm256_mul_epu32(a2, mr2[2]);	// #22
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[3]);	// #23
	ct2 = _mm256_mul_epu32(a2, mr2[3]);	// #24
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	cs1 = _mm256_add_epi64(cs1, cb1);
	cs2 = _mm256_add_epi64(cs2, cb2);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[4]);
	cs2 = _mm256_add_epi64(cs2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #25
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #26
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #27
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #28
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #29
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #30
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #31
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #32
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[3]);	// #33
	cs2 = _mm256_mul_epu32(cm2, mp[3]);	// #34
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[4] = _mm256_add_epi64(c1, t1[5]);
	t2[4] = _mm256_add_epi64(c2, t2[5]);

	a1 = mr1[2];
	a2 = mr2[2];
	cs1 = _mm256_mul_epu32(a1, a1);	// #35
	cs2 = _mm256_mul_epu32(a2, a2);	// #36
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	ct1 = _mm256_mul_epu32(a1, mr1[3]);	// #37
	ct2 = _mm256_mul_epu32(a2, mr2[3]);	// #38
	cs1 = _mm256_add_epi64(ct1, ca1);
	cs2 = _mm256_add_epi64(ct2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	cd1 = _mm256_blend_epi32(cs1, zero, 0xaa);
	cd2 = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ct1, cd1);
	cs2 = _mm256_add_epi64(ct2, cd2);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	cb1 = _mm256_srli_epi64(cs1, 32);
	cb2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(cb1, ca1);
	cs2 = _mm256_add_epi64(cb2, ca2);
	cs1 = _mm256_add_epi64(cs1, t1[4]);
	cs2 = _mm256_add_epi64(cs2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #39
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #40
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #41
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #42
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #43
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #44
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	t1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #45
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #46
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	t1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[3]);	// #47
	cs2 = _mm256_mul_epu32(cm2, mp[3]);	// #48
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	t1[4] = _mm256_add_epi64(c1, t1[5]);
	t2[4] = _mm256_add_epi64(c2, t2[5]);

	a1 = mr1[3];
	a2 = mr2[3];
	cs1 = _mm256_mul_epu32(a1, a1);	// #49
	cs2 = _mm256_mul_epu32(a2, a2);	// #50
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	t1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	ca1 = _mm256_srli_epi64(cs1, 32);
	ca2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(ca1, t1[4]);
	cs2 = _mm256_add_epi64(ca2, t2[4]);
	t1[4] = _mm256_blend_epi32(cs1, zero, 0xaa);
	t2[4] = _mm256_blend_epi32(cs2, zero, 0xaa);
	t1[5] = _mm256_srli_epi64(cs1, 32);
	t2[5] = _mm256_srli_epi64(cs2, 32);

	cm1 = _mm256_mul_epu32(mmagic, t1[0]);	// #51
	cm2 = _mm256_mul_epu32(mmagic, t2[0]);	// #52
	cs1 = _mm256_mul_epu32(cm1, mp[0]);	// #53
	cs2 = _mm256_mul_epu32(cm2, mp[0]);	// #54
	cs1 = _mm256_add_epi64(cs1, t1[0]);
	cs2 = _mm256_add_epi64(cs2, t2[0]);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[1]);	// #55
	cs2 = _mm256_mul_epu32(cm2, mp[1]);	// #56
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[1]);
	cs2 = _mm256_add_epi64(cs2, t2[1]);
	mr1[0] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[0] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[2]);	// #57
	cs2 = _mm256_mul_epu32(cm2, mp[2]);	// #58
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[2]);
	cs2 = _mm256_add_epi64(cs2, t2[2]);
	mr1[1] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[1] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_mul_epu32(cm1, mp[3]);	// #59
	cs2 = _mm256_mul_epu32(cm2, mp[3]);	// #60
	cs1 = _mm256_add_epi64(cs1, c1);
	cs2 = _mm256_add_epi64(cs2, c2);
	cs1 = _mm256_add_epi64(cs1, t1[3]);
	cs2 = _mm256_add_epi64(cs2, t2[3]);
	mr1[2] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[2] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[4]);
	cs2 = _mm256_add_epi64(c2, t2[4]);
	mr1[3] = _mm256_blend_epi32(cs1, zero, 0xaa);
	mr2[3] = _mm256_blend_epi32(cs2, zero, 0xaa);
	c1 = _mm256_srli_epi64(cs1, 32);
	c2 = _mm256_srli_epi64(cs2, 32);
	cs1 = _mm256_add_epi64(c1, t1[5]);
	cs2 = _mm256_add_epi64(c2, t2[5]);

	avx22_modsubtract54(mr1, cs1, mr2, cs2, mp);
}

// -----------------------------------------------------------------------------------

bool avx22_sprp4(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m256i p[4], r1[4], r2[4], one[4], m1[4], b1[4], b2[4], bbb1[4], bbb2[4];
	uint64_t bit, k;
	uint128_t s;
	__m256i mmagic = _mm256_set1_epi64x(mm);
	b1[0] = _mm256_set_epi64x((uint32_t) (bases[0] >> 0),
				  (uint32_t) (bases[1] >> 0), (uint32_t) (bases[2] >> 0), (uint32_t) (bases[3] >> 0));
	b2[0] = _mm256_set_epi64x((uint32_t) (bases[4] >> 0),
				  (uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm256_set1_epi64x((uint32_t) (on >> 0));
	p[0] = _mm256_set1_epi64x((uint32_t) (v >> 0));
	b1[1] = _mm256_set_epi64x((uint32_t) (bases[0] >> 32),
				  (uint32_t) (bases[1] >> 32), (uint32_t) (bases[2] >> 32), (uint32_t) (bases[3] >> 32));
	b2[1] = _mm256_set_epi64x((uint32_t) (bases[4] >> 32),
				  (uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm256_set1_epi64x((uint32_t) (on >> 32));
	p[1] = _mm256_set1_epi64x((uint32_t) (v >> 32));
	b1[2] = _mm256_set_epi64x((uint32_t) (bases[0] >> 64),
				  (uint32_t) (bases[1] >> 64), (uint32_t) (bases[2] >> 64), (uint32_t) (bases[3] >> 64));
	b2[2] = _mm256_set_epi64x((uint32_t) (bases[4] >> 64),
				  (uint32_t) (bases[5] >> 64), (uint32_t) (bases[6] >> 64), (uint32_t) (bases[7] >> 64));
	one[2] = _mm256_set1_epi64x((uint32_t) (on >> 64));
	p[2] = _mm256_set1_epi64x((uint32_t) (v >> 64));
	b1[3] = _mm256_set_epi64x((uint32_t) (bases[0] >> 96),
				  (uint32_t) (bases[1] >> 96), (uint32_t) (bases[2] >> 96), (uint32_t) (bases[3] >> 96));
	b2[3] = _mm256_set_epi64x((uint32_t) (bases[4] >> 96),
				  (uint32_t) (bases[5] >> 96), (uint32_t) (bases[6] >> 96), (uint32_t) (bases[7] >> 96));
	one[3] = _mm256_set1_epi64x((uint32_t) (on >> 96));
	p[3] = _mm256_set1_epi64x((uint32_t) (v >> 96));
// p - 1
	avx22_subtract4(m1, p, one);
// first value
	r1[0] = b1[0];
	r2[0] = b2[0];
	r1[1] = b1[1];
	r2[1] = b2[1];
	r1[2] = b1[2];
	r2[2] = b2[2];
	r1[3] = b1[3];
	r2[3] = b2[3];
// Windowing
	bbb1[0] = b1[0];
	bbb2[0] = b2[0];
	bbb1[1] = b1[1];
	bbb2[1] = b2[1];
	bbb1[2] = b1[2];
	bbb2[2] = b2[2];
	bbb1[3] = b1[3];
	bbb2[3] = b2[3];
	avx22_modsqu4(bbb1, bbb2, p, mmagic);
	avx22_modmul4(bbb1, b1, bbb2, b2, p, mmagic);
// MR exponentiation bit per bit
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx22_modsqu4(r1, r2, p, mmagic);
			avx22_modsqu4(r1, r2, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx22_modsqu4(r1, r2, p, mmagic);
			break;
		case 2:
			avx22_modsqu4(r1, r2, p, mmagic);
			avx22_modmul4(r1, b1, r2, b2, p, mmagic);
			avx22_modsqu4(r1, r2, p, mmagic);
			break;
		case 3:
			avx22_modsqu4(r1, r2, p, mmagic);
			avx22_modsqu4(r1, r2, p, mmagic);
			avx22_modmul4(r1, bbb1, r2, bbb2, p, mmagic);
			break;
		}
	}
	while (bit > 0) {
		bit--;
		// Square
		avx22_modsqu4(r1, r2, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx22_modmul4(r1, b1, r2, b2, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m256i zero[4];
	__m256i mask1 = _mm256_setzero_si256();
	__m256i mask2 = _mm256_setzero_si256();
	zero[0] = _mm256_setzero_si256();
	zero[1] = _mm256_setzero_si256();
	zero[2] = _mm256_setzero_si256();
	zero[3] = _mm256_setzero_si256();
	if (avx22_cmp_next4(&mask1, b1, zero)) {
		return true;
	}
	if (avx22_cmp_next4(&mask2, b2, zero)) {
		return true;
	}
// check current result == 1
	if (avx22_cmp_next4(&mask1, r1, one)) {
		return true;
	}
	if (avx22_cmp_next4(&mask2, r2, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx22_cmp_next4(&mask1, r1, m1)) {
			return true;
		}
		if (avx22_cmp_next4(&mask2, r2, m1)) {
			return true;
		}
		// square
		avx22_modsqu4(r1, r2, p, mmagic);
		// check current result == 1
		if (avx22_neg_cmp_next4(&mask1, r1, one)) {
			// non-trivial quadratic residue found
			return false;
		}
		if (avx22_neg_cmp_next4(&mask2, r2, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx22_cmp_next4(&mask1, r1, m1)) {
		return true;
	}
	if (avx22_cmp_next4(&mask2, r2, m1)) {
		return true;
	}
	return false;
}

#endif				// avx2

// -----------------------------------------------------------------------------------
//
//Generated code ends here
//
// -----------------------------------------------------------------------------------

#if defined(__AVX512F__)

// -----------------------------------------------------------------------------------
//
//Generated code starts here
//
// -----------------------------------------------------------------------------------

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract1(__m512i * mr, __m512i mp[1])
{
	__m512i t[1], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract21(__m512i * mr, __m512i mc, __m512i mp[1])
{
	__m512i t[1], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mc);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_subtract1(__m512i * mr, __m512i ma[1], __m512i mb[1])
{
	mr[0] = _mm512_sub_epi64(ma[0], mb[0]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_cmp_next1(uint64_t * mask, __m512i ma[1], __m512i mb[1])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 |= *mask;
	*mask = t0;
	return t0 == (unsigned)0xff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_neg_cmp_next1(uint64_t * mask, __m512i ma[1], __m512i mb[1])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	uint64_t f = (~*mask) & t0;
	return f != 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modmul1(__m512i * mr, __m512i mb[1], __m512i mp[1], __m512i mmagic)
{
	__m512i t[3];
	__m512i m, c, cs, d, ds, a;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #2
	t[1] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[0]);	// #3
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	ds = _mm512_srli_epi64(ds, 32);
	avx512_modsubtract21(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modsqu1(__m512i * mr, __m512i mp[1], __m512i mmagic)
{
	__m512i t[3];
	__m512i a, m, c, c1, c2, cs, ct, cd;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, a);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #2
	t[1] = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[0]);	// #3
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_srli_epi64(cs, 32);
	avx512_modsubtract21(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx512_sprp1(uint32_t v, uint32_t mm, uint32_t on, uint32_t * bases)
{
	__m512i p[1], r[1], one[1], m1[1], b[1], bbb[1];
	uint64_t bit, k;
	uint32_t s;
	__m512i mmagic = _mm512_set1_epi64(mm);
	b[0] = _mm512_set_epi64((uint32_t) (bases[0] >> 0),
				(uint32_t) (bases[1] >> 0),
				(uint32_t) (bases[2] >> 0),
				(uint32_t) (bases[3] >> 0),
				(uint32_t) (bases[4] >> 0),
				(uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm512_set1_epi64((uint32_t) (on >> 0));
	p[0] = _mm512_set1_epi64((uint32_t) (v >> 0));
// p - 1
	avx512_subtract1(m1, p, one);
// window bbb = b^3 mod p
	bbb[0] = b[0];
	avx512_modsqu1(bbb, p, mmagic);
	avx512_modmul1(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
// exponentiation bit per bit with sliding window size 2
	k = my_ctz32(v - 1);
	s = v >> k;
	bit = 31 - my_clz32(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx512_modsqu1(r, p, mmagic);
			avx512_modsqu1(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx512_modsqu1(r, p, mmagic);
			break;
		case 2:
			avx512_modsqu1(r, p, mmagic);
			avx512_modmul1(r, b, p, mmagic);
			avx512_modsqu1(r, p, mmagic);
			break;
		case 3:
			avx512_modsqu1(r, p, mmagic);
			avx512_modsqu1(r, p, mmagic);
			avx512_modmul1(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
	while (bit > 0) {
		bit--;
		// Square
		avx512_modsqu1(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx512_modmul1(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m512i zero[1];
	uint64_t mask = 0;
	zero[0] = _mm512_setzero_si512();
	if (avx512_cmp_next1(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx512_cmp_next1(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx512_cmp_next1(&mask, r, m1)) {
			return true;
		}
		// square
		avx512_modsqu1(r, p, mmagic);
		// check current result == 1
		if (avx512_neg_cmp_next1(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx512_cmp_next1(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract2(__m512i * mr, __m512i mp[2])
{
	__m512i t[2], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract32(__m512i * mr, __m512i mc, __m512i mp[2])
{
	__m512i t[2], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mc);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_subtract2(__m512i * mr, __m512i ma[2], __m512i mb[2])
{
	__m512i cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[1]);
	mr[1] = _mm512_sub_epi64(cs, mb[1]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_cmp_next2(uint64_t * mask, __m512i ma[2], __m512i mb[2])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	t0 |= *mask;
	*mask = t0;
	return t0 == (unsigned)0xff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_neg_cmp_next2(uint64_t * mask, __m512i ma[2], __m512i mb[2])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	uint64_t f = (~*mask) & t0;
	return f != 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modmul2(__m512i * mr, __m512i mb[2], __m512i mp[2], __m512i mmagic)
{
	__m512i t[4];
	__m512i m, c, cs, d, ds, a;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm512_mul_epu32(m, mp[0]);	// #3
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #4
	cs = _mm512_add_epi64(cs, c);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#5
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, c);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	t[2] = _mm512_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm512_mul_epu32(a, mb[0]);	// #6
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #7
	ds = _mm512_mul_epu32(m, mp[0]);	// #8
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #9
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#10
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[3] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, t[3]);

	avx512_modsubtract32(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modsqu2(__m512i * mr, __m512i mp[2], __m512i mmagic)
{
	__m512i t[4];
	__m512i a, m, c, c1, c2, cs, ct, cd;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, a);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[1]);	// #2
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[3] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #3
	cs = _mm512_mul_epu32(m, mp[0]);	// #4
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #5
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[2] = _mm512_add_epi64(c, t[3]);

	a = mr[1];
	cs = _mm512_mul_epu32(a, a);	// #6
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c1, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[3] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #7
	cs = _mm512_mul_epu32(m, mp[0]);	// #8
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #9
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[3]);

	avx512_modsubtract32(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx512_sprp2(uint64_t v, uint32_t mm, uint64_t on, uint64_t * bases)
{
	__m512i p[2], r[2], one[2], m1[2], b[2], bbb[2];
	uint64_t bit, k;
	uint64_t s;
	__m512i mmagic = _mm512_set1_epi64(mm);
	b[0] = _mm512_set_epi64((uint32_t) (bases[0] >> 0),
				(uint32_t) (bases[1] >> 0),
				(uint32_t) (bases[2] >> 0),
				(uint32_t) (bases[3] >> 0),
				(uint32_t) (bases[4] >> 0),
				(uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm512_set1_epi64((uint32_t) (on >> 0));
	p[0] = _mm512_set1_epi64((uint32_t) (v >> 0));
	b[1] = _mm512_set_epi64((uint32_t) (bases[0] >> 32),
				(uint32_t) (bases[1] >> 32),
				(uint32_t) (bases[2] >> 32),
				(uint32_t) (bases[3] >> 32),
				(uint32_t) (bases[4] >> 32),
				(uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm512_set1_epi64((uint32_t) (on >> 32));
	p[1] = _mm512_set1_epi64((uint32_t) (v >> 32));
// p - 1
	avx512_subtract2(m1, p, one);
// window bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	avx512_modsqu2(bbb, p, mmagic);
	avx512_modmul2(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
// exponentiation bit per bit with sliding window size 2
	k = my_ctz64(v - 1);
	s = v >> k;
	bit = 63 - my_clz64(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx512_modsqu2(r, p, mmagic);
			avx512_modsqu2(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx512_modsqu2(r, p, mmagic);
			break;
		case 2:
			avx512_modsqu2(r, p, mmagic);
			avx512_modmul2(r, b, p, mmagic);
			avx512_modsqu2(r, p, mmagic);
			break;
		case 3:
			avx512_modsqu2(r, p, mmagic);
			avx512_modsqu2(r, p, mmagic);
			avx512_modmul2(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
	while (bit > 0) {
		bit--;
		// Square
		avx512_modsqu2(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx512_modmul2(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m512i zero[2];
	uint64_t mask = 0;
	zero[0] = _mm512_setzero_si512();
	zero[1] = _mm512_setzero_si512();
	if (avx512_cmp_next2(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx512_cmp_next2(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx512_cmp_next2(&mask, r, m1)) {
			return true;
		}
		// square
		avx512_modsqu2(r, p, mmagic);
		// check current result == 1
		if (avx512_neg_cmp_next2(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx512_cmp_next2(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract3(__m512i * mr, __m512i mp[3])
{
	__m512i t[3], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[2]);
	cs = _mm512_sub_epi64(cs, mp[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
	mr[2] = _mm512_mask_shuffle_epi32(t[2], f, mr[2], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract43(__m512i * mr, __m512i mc, __m512i mp[3])
{
	__m512i t[3], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[2]);
	cs = _mm512_sub_epi64(cs, mp[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mc);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
	mr[2] = _mm512_mask_shuffle_epi32(t[2], f, mr[2], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_subtract3(__m512i * mr, __m512i ma[3], __m512i mb[3])
{
	__m512i cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[1]);
	cs = _mm512_sub_epi64(cs, mb[1]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[2]);
	mr[2] = _mm512_sub_epi64(cs, mb[2]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_cmp_next3(uint64_t * mask, __m512i ma[3], __m512i mb[3])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[2], mb[2]));
	t0 |= *mask;
	*mask = t0;
	return t0 == (unsigned)0xff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_neg_cmp_next3(uint64_t * mask, __m512i ma[3], __m512i mb[3])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[2], mb[2]));
	uint64_t f = (~*mask) & t0;
	return f != 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modmul3(__m512i * mr, __m512i mb[3], __m512i mp[3], __m512i mmagic)
{
	__m512i t[5];
	__m512i m, c, cs, d, ds, a;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm512_mul_epu32(m, mp[0]);	// #3
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #4
	cs = _mm512_add_epi64(cs, c);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#5
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #6
	cs = _mm512_add_epi64(cs, c);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#7
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, c);
	t[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	t[3] = _mm512_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm512_mul_epu32(a, mb[0]);	// #8
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #9
	ds = _mm512_mul_epu32(m, mp[0]);	// #10
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #11
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#12
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #13
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#14
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	t[3] = _mm512_add_epi64(d, t[4]);

	a = mr[2];
	cs = _mm512_mul_epu32(a, mb[0]);	// #15
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #16
	ds = _mm512_mul_epu32(m, mp[0]);	// #17
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #18
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#19
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #20
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#21
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[3]);
	mr[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, t[4]);

	avx512_modsubtract43(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modsqu3(__m512i * mr, __m512i mp[3], __m512i mmagic)
{
	__m512i t[5];
	__m512i a, m, c, c1, c2, cs, ct, cd;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, a);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[1]);	// #2
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[2]);	// #3
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	cs = _mm512_add_epi64(cs, c2);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #4
	cs = _mm512_mul_epu32(m, mp[0]);	// #5
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #6
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #7
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[3] = _mm512_add_epi64(c, t[4]);

	a = mr[1];
	cs = _mm512_mul_epu32(a, a);	// #8
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[2]);	// #9
	cs = _mm512_add_epi64(ct, c1);
	cs = _mm512_add_epi64(cs, t[2]);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	cs = _mm512_add_epi64(cs, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #10
	cs = _mm512_mul_epu32(m, mp[0]);	// #11
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #12
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #13
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[3] = _mm512_add_epi64(c, t[4]);

	a = mr[2];
	cs = _mm512_mul_epu32(a, a);	// #14
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c1, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #15
	cs = _mm512_mul_epu32(m, mp[0]);	// #16
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #17
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #18
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[3]);
	mr[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[4]);

	avx512_modsubtract43(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx512_sprp3(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m512i p[3], r[3], one[3], m1[3], b[3], bbb[3];
	uint64_t bit, k;
	uint128_t s;
	__m512i mmagic = _mm512_set1_epi64(mm);
	b[0] = _mm512_set_epi64((uint32_t) (bases[0] >> 0),
				(uint32_t) (bases[1] >> 0),
				(uint32_t) (bases[2] >> 0),
				(uint32_t) (bases[3] >> 0),
				(uint32_t) (bases[4] >> 0),
				(uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm512_set1_epi64((uint32_t) (on >> 0));
	p[0] = _mm512_set1_epi64((uint32_t) (v >> 0));
	b[1] = _mm512_set_epi64((uint32_t) (bases[0] >> 32),
				(uint32_t) (bases[1] >> 32),
				(uint32_t) (bases[2] >> 32),
				(uint32_t) (bases[3] >> 32),
				(uint32_t) (bases[4] >> 32),
				(uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm512_set1_epi64((uint32_t) (on >> 32));
	p[1] = _mm512_set1_epi64((uint32_t) (v >> 32));
	b[2] = _mm512_set_epi64((uint32_t) (bases[0] >> 64),
				(uint32_t) (bases[1] >> 64),
				(uint32_t) (bases[2] >> 64),
				(uint32_t) (bases[3] >> 64),
				(uint32_t) (bases[4] >> 64),
				(uint32_t) (bases[5] >> 64), (uint32_t) (bases[6] >> 64), (uint32_t) (bases[7] >> 64));
	one[2] = _mm512_set1_epi64((uint32_t) (on >> 64));
	p[2] = _mm512_set1_epi64((uint32_t) (v >> 64));
// p - 1
	avx512_subtract3(m1, p, one);
// window bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	bbb[2] = b[2];
	avx512_modsqu3(bbb, p, mmagic);
	avx512_modmul3(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
// exponentiation bit per bit with sliding window size 2
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx512_modsqu3(r, p, mmagic);
			avx512_modsqu3(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx512_modsqu3(r, p, mmagic);
			break;
		case 2:
			avx512_modsqu3(r, p, mmagic);
			avx512_modmul3(r, b, p, mmagic);
			avx512_modsqu3(r, p, mmagic);
			break;
		case 3:
			avx512_modsqu3(r, p, mmagic);
			avx512_modsqu3(r, p, mmagic);
			avx512_modmul3(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
	while (bit > 0) {
		bit--;
		// Square
		avx512_modsqu3(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx512_modmul3(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m512i zero[3];
	uint64_t mask = 0;
	zero[0] = _mm512_setzero_si512();
	zero[1] = _mm512_setzero_si512();
	zero[2] = _mm512_setzero_si512();
	if (avx512_cmp_next3(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx512_cmp_next3(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx512_cmp_next3(&mask, r, m1)) {
			return true;
		}
		// square
		avx512_modsqu3(r, p, mmagic);
		// check current result == 1
		if (avx512_neg_cmp_next3(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx512_cmp_next3(&mask, r, m1)) {
		return true;
	}
	return false;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract4(__m512i * mr, __m512i mp[4])
{
	__m512i t[4], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[2]);
	cs = _mm512_sub_epi64(cs, mp[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[3]);
	cs = _mm512_sub_epi64(cs, mp[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
	mr[2] = _mm512_mask_shuffle_epi32(t[2], f, mr[2], _MM_PERM_DCBA);
	mr[3] = _mm512_mask_shuffle_epi32(t[3], f, mr[3], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_modsubtract54(__m512i * mr, __m512i mc, __m512i mp[4])
{
	__m512i t[4], cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(mr[0], mp[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[1]);
	cs = _mm512_sub_epi64(cs, mp[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[2]);
	cs = _mm512_sub_epi64(cs, mp[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mr[3]);
	cs = _mm512_sub_epi64(cs, mp[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, mc);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	__mmask16 f = _mm512_movepi32_mask(cs);
	mr[0] = _mm512_mask_shuffle_epi32(t[0], f, mr[0], _MM_PERM_DCBA);
	mr[1] = _mm512_mask_shuffle_epi32(t[1], f, mr[1], _MM_PERM_DCBA);
	mr[2] = _mm512_mask_shuffle_epi32(t[2], f, mr[2], _MM_PERM_DCBA);
	mr[3] = _mm512_mask_shuffle_epi32(t[3], f, mr[3], _MM_PERM_DCBA);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
void avx512_subtract4(__m512i * mr, __m512i ma[4], __m512i mb[4])
{
	__m512i cs;
	__mmask16 k = _mm512_int2mask(0x5555);
	cs = _mm512_sub_epi64(ma[0], mb[0]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[1]);
	cs = _mm512_sub_epi64(cs, mb[1]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[2]);
	cs = _mm512_sub_epi64(cs, mb[2]);
	mr[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	cs = _mm512_shuffle_epi32(cs, _MM_PERM_DDBB);
	cs = _mm512_add_epi64(cs, ma[3]);
	mr[3] = _mm512_sub_epi64(cs, mb[3]);
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_cmp_next4(uint64_t * mask, __m512i ma[4], __m512i mb[4])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[2], mb[2]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[3], mb[3]));
	t0 |= *mask;
	*mask = t0;
	return t0 == (unsigned)0xff;
}

// -----------------------------------------------------------------------------------

static inline __attribute__((always_inline))
bool avx512_neg_cmp_next4(uint64_t * mask, __m512i ma[4], __m512i mb[4])
{
	uint64_t t0;
	t0 = _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[0], mb[0]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[1], mb[1]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[2], mb[2]));
	t0 &= _mm512_mask2int(_mm512_cmpeq_epu64_mask(ma[3], mb[3]));
	uint64_t f = (~*mask) & t0;
	return f != 0x0;
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modmul4(__m512i * mr, __m512i mb[4], __m512i mp[4], __m512i mmagic)
{
	__m512i t[6];
	__m512i m, c, cs, d, ds, a;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, mb[0]);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #2
	ds = _mm512_mul_epu32(m, mp[0]);	// #3
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #4
	cs = _mm512_add_epi64(cs, c);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#5
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #6
	cs = _mm512_add_epi64(cs, c);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#7
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[3]);	// #8
	cs = _mm512_add_epi64(cs, c);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[3]);	//#9
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, c);
	t[3] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	t[4] = _mm512_srli_epi64(ds, 32);

	a = mr[1];
	cs = _mm512_mul_epu32(a, mb[0]);	// #10
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #11
	ds = _mm512_mul_epu32(m, mp[0]);	// #12
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #13
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#14
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #15
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#16
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[3]);	// #17
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[3]);	//#18
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[4]);
	t[3] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	t[4] = _mm512_add_epi64(d, t[5]);

	a = mr[2];
	cs = _mm512_mul_epu32(a, mb[0]);	// #19
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #20
	ds = _mm512_mul_epu32(m, mp[0]);	// #21
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #22
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#23
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #24
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#25
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[3]);	// #26
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[3]);	//#27
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[4]);
	t[3] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	t[4] = _mm512_add_epi64(d, t[5]);

	a = mr[3];
	cs = _mm512_mul_epu32(a, mb[0]);	// #28
	cs = _mm512_add_epi64(cs, t[0]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	m = _mm512_mul_epu32(mmagic, t[0]);	// #29
	ds = _mm512_mul_epu32(m, mp[0]);	// #30
	ds = _mm512_add_epi64(ds, t[0]);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[1]);	// #31
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[1]);	//#32
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[2]);	// #33
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[2]);	//#34
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_mul_epu32(a, mb[3]);	// #35
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	ds = _mm512_mul_epu32(m, mp[3]);	//#36
	ds = _mm512_add_epi64(ds, d);
	ds = _mm512_add_epi64(ds, t[3]);
	mr[2] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);
	ds = _mm512_add_epi64(d, t[4]);
	mr[3] = _mm512_maskz_shuffle_epi32(k, ds, _MM_PERM_DCBA);
	d = _mm512_srli_epi64(ds, 32);
	ds = _mm512_add_epi64(d, t[5]);

	avx512_modsubtract54(mr, ds, mp);
}

// -----------------------------------------------------------------------------------

static
    inline __attribute__((always_inline))
void avx512_modsqu4(__m512i * mr, __m512i mp[4], __m512i mmagic)
{
	__m512i t[6];
	__m512i a, m, c, c1, c2, cs, ct, cd;
	__mmask16 k = 0x5555;

	a = mr[0];
	cs = _mm512_mul_epu32(a, a);	// #1
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[1]);	// #2
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[2]);	// #3
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	cs = _mm512_add_epi64(cs, c2);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[3]);	// #4
	cs = _mm512_add_epi64(ct, c1);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	cs = _mm512_add_epi64(cs, c2);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #5
	cs = _mm512_mul_epu32(m, mp[0]);	// #6
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #7
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #8
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[3]);	// #9
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[4] = _mm512_add_epi64(c, t[5]);

	a = mr[1];
	cs = _mm512_mul_epu32(a, a);	// #10
	cs = _mm512_add_epi64(cs, t[1]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[2]);	// #11
	cs = _mm512_add_epi64(ct, c1);
	cs = _mm512_add_epi64(cs, t[2]);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[3]);	// #12
	cs = _mm512_add_epi64(ct, c1);
	cs = _mm512_add_epi64(cs, t[3]);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	cs = _mm512_add_epi64(cs, c2);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	cs = _mm512_add_epi64(cs, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #13
	cs = _mm512_mul_epu32(m, mp[0]);	// #14
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #15
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #16
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[3]);	// #17
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[4] = _mm512_add_epi64(c, t[5]);

	a = mr[2];
	cs = _mm512_mul_epu32(a, a);	// #18
	cs = _mm512_add_epi64(cs, t[2]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	ct = _mm512_mul_epu32(a, mr[3]);	// #19
	cs = _mm512_add_epi64(ct, c1);
	cs = _mm512_add_epi64(cs, t[3]);
	cd = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(ct, cd);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c2 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c2, c1);
	cs = _mm512_add_epi64(cs, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #20
	cs = _mm512_mul_epu32(m, mp[0]);	// #21
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #22
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	t[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #23
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	t[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[3]);	// #24
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	t[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[4]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	t[4] = _mm512_add_epi64(c, t[5]);

	a = mr[3];
	cs = _mm512_mul_epu32(a, a);	// #25
	cs = _mm512_add_epi64(cs, t[3]);
	t[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c1 = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c1, t[4]);
	t[4] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	t[5] = _mm512_srli_epi64(cs, 32);

	m = _mm512_mul_epu32(mmagic, t[0]);	// #26
	cs = _mm512_mul_epu32(m, mp[0]);	// #27
	cs = _mm512_add_epi64(cs, t[0]);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[1]);	// #28
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[1]);
	mr[0] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[2]);	// #29
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[2]);
	mr[1] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_mul_epu32(m, mp[3]);	// #30
	cs = _mm512_add_epi64(cs, c);
	cs = _mm512_add_epi64(cs, t[3]);
	mr[2] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[4]);
	mr[3] = _mm512_maskz_shuffle_epi32(k, cs, _MM_PERM_DCBA);
	c = _mm512_srli_epi64(cs, 32);
	cs = _mm512_add_epi64(c, t[5]);

	avx512_modsubtract54(mr, cs, mp);
}

// -----------------------------------------------------------------------------------

bool avx512_sprp4(uint128_t v, uint32_t mm, uint128_t on, uint128_t * bases)
{
	__m512i p[4], r[4], one[4], m1[4], b[4], bbb[4];
	uint64_t bit, k;
	uint128_t s;
	__m512i mmagic = _mm512_set1_epi64(mm);
	b[0] = _mm512_set_epi64((uint32_t) (bases[0] >> 0),
				(uint32_t) (bases[1] >> 0),
				(uint32_t) (bases[2] >> 0),
				(uint32_t) (bases[3] >> 0),
				(uint32_t) (bases[4] >> 0),
				(uint32_t) (bases[5] >> 0), (uint32_t) (bases[6] >> 0), (uint32_t) (bases[7] >> 0));
	one[0] = _mm512_set1_epi64((uint32_t) (on >> 0));
	p[0] = _mm512_set1_epi64((uint32_t) (v >> 0));
	b[1] = _mm512_set_epi64((uint32_t) (bases[0] >> 32),
				(uint32_t) (bases[1] >> 32),
				(uint32_t) (bases[2] >> 32),
				(uint32_t) (bases[3] >> 32),
				(uint32_t) (bases[4] >> 32),
				(uint32_t) (bases[5] >> 32), (uint32_t) (bases[6] >> 32), (uint32_t) (bases[7] >> 32));
	one[1] = _mm512_set1_epi64((uint32_t) (on >> 32));
	p[1] = _mm512_set1_epi64((uint32_t) (v >> 32));
	b[2] = _mm512_set_epi64((uint32_t) (bases[0] >> 64),
				(uint32_t) (bases[1] >> 64),
				(uint32_t) (bases[2] >> 64),
				(uint32_t) (bases[3] >> 64),
				(uint32_t) (bases[4] >> 64),
				(uint32_t) (bases[5] >> 64), (uint32_t) (bases[6] >> 64), (uint32_t) (bases[7] >> 64));
	one[2] = _mm512_set1_epi64((uint32_t) (on >> 64));
	p[2] = _mm512_set1_epi64((uint32_t) (v >> 64));
	b[3] = _mm512_set_epi64((uint32_t) (bases[0] >> 96),
				(uint32_t) (bases[1] >> 96),
				(uint32_t) (bases[2] >> 96),
				(uint32_t) (bases[3] >> 96),
				(uint32_t) (bases[4] >> 96),
				(uint32_t) (bases[5] >> 96), (uint32_t) (bases[6] >> 96), (uint32_t) (bases[7] >> 96));
	one[3] = _mm512_set1_epi64((uint32_t) (on >> 96));
	p[3] = _mm512_set1_epi64((uint32_t) (v >> 96));
// p - 1
	avx512_subtract4(m1, p, one);
// window bbb = b^3 mod p
	bbb[0] = b[0];
	bbb[1] = b[1];
	bbb[2] = b[2];
	bbb[3] = b[3];
	avx512_modsqu4(bbb, p, mmagic);
	avx512_modmul4(bbb, b, p, mmagic);
// first value
	r[0] = b[0];
	r[1] = b[1];
	r[2] = b[2];
	r[3] = b[3];
// exponentiation bit per bit with sliding window size 2
	k = my_ctz128(v - 1);
	s = v >> k;
	bit = 127 - my_clz128(s);
	while (bit > 1) {
		bit -= 2;
		switch ((s >> bit) & 3) {
		case 0:
			avx512_modsqu4(r, p, mmagic);
			avx512_modsqu4(r, p, mmagic);
			break;
		case 1:
			bit += 1;
			avx512_modsqu4(r, p, mmagic);
			break;
		case 2:
			avx512_modsqu4(r, p, mmagic);
			avx512_modmul4(r, b, p, mmagic);
			avx512_modsqu4(r, p, mmagic);
			break;
		case 3:
			avx512_modsqu4(r, p, mmagic);
			avx512_modsqu4(r, p, mmagic);
			avx512_modmul4(r, bbb, p, mmagic);
			break;
		}
	}
// exponentiation bit per bit
	while (bit > 0) {
		bit--;
		// Square
		avx512_modsqu4(r, p, mmagic);
		if ((s >> bit) & 1) {
			// Multiply
			avx512_modmul4(r, b, p, mmagic);
		}
	}
// check bases which are 0 mod n (they must return true)
	__m512i zero[4];
	uint64_t mask = 0;
	zero[0] = _mm512_setzero_si512();
	zero[1] = _mm512_setzero_si512();
	zero[2] = _mm512_setzero_si512();
	zero[3] = _mm512_setzero_si512();
	if (avx512_cmp_next4(&mask, b, zero)) {
		return true;
	}
// check current result == 1
	if (avx512_cmp_next4(&mask, r, one)) {
		return true;
	}
// MR iteration square per square
	while (k > 1) {
		k -= 1;
		// check current result == m-1
		if (avx512_cmp_next4(&mask, r, m1)) {
			return true;
		}
		// square
		avx512_modsqu4(r, p, mmagic);
		// check current result == 1
		if (avx512_neg_cmp_next4(&mask, r, one)) {
			// non-trivial quadratic residue found
			return false;
		}
	}
// check current result == m-1
	if (avx512_cmp_next4(&mask, r, m1)) {
		return true;
	}
	return false;
}

#endif				// avx2

// -----------------------------------------------------------------------------------
//
//Generated code ends here
//
// -----------------------------------------------------------------------------------

// precomputations for montgomery implementation based on 32x32 avx multiplication instructions (mul_epu32)

#define BASES(type, count, mod, var, base)\
	type var = base;  \
	do { while (var >= mod) var -= mod;   \
	montg_bases[count] = var; \
	} while (0)

#define BASES3(type, count, mod, var, base1, base2, base3)\
	type var = base1;  \
	do { var += base2;  \
	if (var < base2) var -= mod; \
	var += base3;  \
	if (var < base3) var -= mod; \
	while (var >= mod) var -= mod;   \
	montg_bases[count] = var; \
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

	BASES(uint64_t, 4, v, t13, t11 + t2);
	BASES(uint64_t, 5, v, t17, t11 + t5 + t1);
	BASES(uint64_t, 6, v, t19, t17 + t2);
	BASES(uint64_t, 7, v, t23, t11 + t11 + t1);
	if (count <= 8)
		return t1;

	BASES(uint64_t, 8, v, t29, t11 + t17 + t1);
	BASES(uint64_t, 9, v, t31, t29 + t2);
	BASES(uint64_t, 10, v, t37, t17 + t17 + t3);
	BASES(uint64_t, 11, v, t41, t17 + t17 + t7);
	if (count <= 12)
		return t1;

	return t1;
}

uint64_t montgomery_bases2(uint64_t * montg_bases, uint64_t v, uint64_t count)
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

	BASES(uint128_t, 12, v, t43, t41 + t2);
	BASES(uint128_t, 13, v, t47, t23 + t23 + t1);
	BASES(uint128_t, 14, v, t53, t23 + t23 + t7);
	BASES(uint128_t, 15, v, t59, t23 + t23 + t13);
	if (count <= 16)
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

	BASES(uint128_t, 20, v, t79, t29 + t31 + t19);
	BASES(uint128_t, 21, v, t83, t29 + t31 + t23);
	BASES(uint128_t, 22, v, t89, t29 + t31 + t29);
	BASES(uint128_t, 23, v, t97, t29 + t31 + t37);
	if (count <= 24)
		return t1;

	return t1;
}

uint128_t montgomery_bases4(uint128_t * montg_bases, uint128_t v, uint64_t count)
{
	// compute x * 2^128 for different values of x

	uint128_t t1 = -v;
	t1 %= v;		// t1

	BASES3(uint128_t, 0, v, t3, t1, t1, t1);
	BASES3(uint128_t, 1, v, t5, t1, t1, t1);
	BASES3(uint128_t, 2, v, t7, t5, t1, t1);
	BASES3(uint128_t, 3, v, t11, t5, t5, t1);
	if (count <= 4)
		return t1;

	BASES3(uint128_t, 4, v, t13, t11, t1, t1);
	BASES3(uint128_t, 5, v, t17, t11, t5, t1);
	BASES3(uint128_t, 6, v, t19, t17, t1, t1);
	BASES3(uint128_t, 7, v, t23, t11, t11, t1);
	if (count <= 8)
		return t1;

	BASES3(uint128_t, 8, v, t29, t11, t17, t1);
	BASES3(uint128_t, 9, v, t31, t29, t1, t1);
	BASES3(uint128_t, 10, v, t37, t17, t17, t3);
	BASES3(uint128_t, 11, v, t41, t17, t17, t7);
	if (count <= 12)
		return t1;

	BASES3(uint128_t, 12, v, t43, t41, t1, t1);
	BASES3(uint128_t, 13, v, t47, t23, t23, t1);
	BASES3(uint128_t, 14, v, t53, t23, t23, t7);
	BASES3(uint128_t, 15, v, t59, t23, t23, t13);
	if (count <= 16)
		return t1;

	BASES3(uint128_t, 16, v, t61, t29, t31, t1);
	BASES3(uint128_t, 17, v, t67, t29, t31, t7);
	BASES3(uint128_t, 18, v, t71, t29, t31, t11);
	BASES3(uint128_t, 19, v, t73, t29, t31, t13);
	if (count <= 20)
		return t1;

	BASES3(uint128_t, 20, v, t79, t29, t31, t19);
	BASES3(uint128_t, 21, v, t83, t29, t31, t23);
	BASES3(uint128_t, 22, v, t89, t29, t31, t29);
	BASES3(uint128_t, 23, v, t97, t29, t31, t37);
	if (count <= 24)
		return t1;

	return t1;
}

unsigned deterministicSprpCount(uint128_t v)
{
	if (v >> 32 == 0) {
		return 4;
	}
	if (v >> 64 == 0) {
		if (v < 2152302898747ull) {
			return 4;
		}
		if (v < 3825123056546413051ull) {
			return 8;
		}
		return 12;
	}
	if (v >> 96 == 0) {
		// if (v < 3317044064679887385961981)
		if (v < ((uint128_t) 0x2BE69ull << 64) + 0x51ADC5B22410A5FDull) {
			return 12;
		}
		return 16;
	}
	if (v >> 112 == 0) {
		return 16;
	}
	return 20;
}

bool avxSprpTest(uint64_t v_lo, uint64_t v_hi)
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

	unsigned count = deterministicSprpCount(v);

	if (v >> 64 == 0) {
		uint64_t montg_bases[12];
		uint64_t one = montgomery_bases2(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 8) {
#ifdef __AVX512F__
			if (!avx512_sprp2(v, m, one, &montg_bases[i]))
				return false;
#else
			if (!avx22_sprp2(v, m, one, &montg_bases[i]))
				return false;
#endif
			i += 8;
			count -= 8;
		}
		while (count >= 4) {
			if (!avx2_sprp2(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	}
	if (v >> 96 == 0) {
		uint128_t montg_bases[16];
		uint128_t one = montgomery_bases3(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 8) {
#ifdef __AVX512F__
			if (!avx512_sprp3(v, m, one, &montg_bases[i]))
				return false;
#else
			if (!avx22_sprp3(v, m, one, &montg_bases[i]))
				return false;
#endif
			i += 8;
			count -= 8;
		}
		while (count >= 4) {
			if (!avx2_sprp3(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	} else {
		uint128_t montg_bases[20];
		uint128_t one = montgomery_bases4(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 8) {
#ifdef __AVX512F__
			if (!avx512_sprp4(v, m, one, &montg_bases[i]))
				return false;
#else
			if (!avx22_sprp4(v, m, one, &montg_bases[i]))
				return false;
#endif
			i += 8;
			count -= 8;
		}
		while (count >= 4) {
			if (!avx2_sprp4(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	}
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

	unsigned count = deterministicSprpCount(v);

	if (v >> 64 == 0) {
		uint64_t montg_bases[12];
		uint64_t one = montgomery_bases2(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 4) {
			if (!avx2_sprp2(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	}
	if (v >> 96 == 0) {
		uint128_t montg_bases[16];
		uint128_t one = montgomery_bases3(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 4) {
			if (!avx2_sprp3(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	} else {
		uint128_t montg_bases[20];
		uint128_t one = montgomery_bases4(montg_bases, v, count);
		unsigned i = 0;

		while (count >= 4) {
			if (!avx2_sprp4(v, m, one, &montg_bases[i]))
				return false;
			i += 4;
			count -= 4;
		}
		return true;
	}
}

