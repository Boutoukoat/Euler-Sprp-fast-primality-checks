
#include <stdint.h>
#include <m_reg.h>

// ------------------------------------------------------------------------
// fast divisibility of 128 bits numbers by small constants 13, 31, 41, 61 .... 
//
// // assuming that divisibility by 2,3,5,11 has already been done
//
// Modular operations with constant modulus have many more tricks.
//
// Gcc and Clang generate too much code to compute the exact modulus, while
// the only thing needed is the exact divisibility.
// if (mod_mersenne_36 % 109 == 0) return false;        // composite for sure
//
// The simplest equivalent, the Barrett way, might be
//
// if ((uint64_t)(n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull) return false; // divisible by 109
//
// ------------------------------------------------------------------------
//
#if 0
if ((uint64_t) (n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
	return false;		// divisible by 3
if ((uint64_t) (n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
	return false;		// divisible by 5
if ((uint64_t) (n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
	return false;		// divisible by 7
if ((uint64_t) (n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
	return false;		// divisible by 11
if ((uint64_t) (n * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
	return false;		// divisible by 13
if ((uint64_t) (n * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
	return false;		// divisible by 17
if ((uint64_t) (n * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
	return false;		// divisible by 19
if ((uint64_t) (n * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
	return false;		// divisible by 23
if ((uint64_t) (n * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
	return false;		// divisible by 29
if ((uint64_t) (n * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
	return false;		// divisible by 31
if ((uint64_t) (n * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
	return false;		// divisible by 37
if ((uint64_t) (n * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
	return false;		// divisible by 41
if ((uint64_t) (n * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
	return false;		// divisible by 43
if ((uint64_t) (n * 0x51b3bea3677d46cfull) <= 0x0572620ae4c415c9ull)
	return false;		// divisible by 47
if ((uint64_t) (n * 0x21cfb2b78c13521dull) <= 0x04d4873ecade304dull)
	return false;		// divisible by 53
if ((uint64_t) (n * 0xcbeea4e1a08ad8f3ull) <= 0x0456c797dd49c341ull)
	return false;		// divisible by 59
if ((uint64_t) (n * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
	return false;		// divisible by 61
if ((uint64_t) (n * 0xf0b7672a07a44c6bull) <= 0x03d226357e16ece5ull)
	return false;		// divisible by 67
if ((uint64_t) (n * 0x193d4bb7e327a977ull) <= 0x039b0ad12073615aull)
	return false;		// divisible by 71
if ((uint64_t) (n * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
	return false;		// divisible by 73
if ((uint64_t) (n * 0x9b8b577e613716afull) <= 0x033d91d2a2067b23ull)
	return false;		// divisible by 79
if ((uint64_t) (n * 0xa3784a062b2e43dbull) <= 0x03159721ed7e7534ull)
	return false;		// divisible by 83
if ((uint64_t) (n * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
	return false;		// divisible by 89
if ((uint64_t) (n * 0xa3a0fd5c5f02a3a1ull) <= 0x02a3a0fd5c5f02a3ull)
	return false;		// divisible by 97
if ((uint64_t) (n * 0x3a4c0a237c32b16dull) <= 0x0288df0cac5b3f5dull)
	return false;		// divisible by 101
if ((uint64_t) (n * 0xdab7ec1dd3431b57ull) <= 0x027c45979c95204full)
	return false;		// divisible by 103
if ((uint64_t) (n * 0x77a04c8f8d28ac43ull) <= 0x02647c69456217ecull)
	return false;		// divisible by 107
if ((uint64_t) (n * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
	return false;		// divisible by 109
if ((uint64_t) (n * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
	return false;		// divisible by 113
if ((uint64_t) (n * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
	return false;		// divisible by 127
if ((uint64_t) (n * 0x03e88cb3c9484e2bull) <= 0x01f44659e4a42715ull)
	return false;		// divisible by 131
if ((uint64_t) (n * 0xe21a291c077975b9ull) <= 0x01de5d6e3f8868a4ull)
	return false;		// divisible by 137
if ((uint64_t) (n * 0x3aef6ca970586723ull) <= 0x01d77b654b82c339ull)
	return false;		// divisible by 139
if ((uint64_t) (n * 0xdf5b0f768ce2cabdull) <= 0x01b7d6c3dda338b2ull)
	return false;		// divisible by 149
if ((uint64_t) (n * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
	return false;		// divisible by 151
if ((uint64_t) (n * 0x5b4fe5e92c0685b5ull) <= 0x01a16d3f97a4b01aull)
	return false;		// divisible by 157
if ((uint64_t) (n * 0x1f693a1c451ab30bull) <= 0x01920fb49d0e228dull)
	return false;		// divisible by 163
if ((uint64_t) (n * 0x8d07aa27db35a717ull) <= 0x01886e5f0abb0499ull)
	return false;		// divisible by 167
if ((uint64_t) (n * 0x882383b30d516325ull) <= 0x017ad2208e0ecc35ull)
	return false;		// divisible by 173
if ((uint64_t) (n * 0xed6866f8d962ae7bull) <= 0x016e1f76b4337c6cull)
	return false;		// divisible by 179
if ((uint64_t) (n * 0x3454dca410f8ed9dull) <= 0x016a13cd15372904ull)
	return false;		// divisible by 181
if ((uint64_t) (n * 0x1d7ca632ee936f3full) <= 0x01571ed3c506b39aull)
	return false;		// divisible by 191
if ((uint64_t) (n * 0x70bf015390948f41ull) <= 0x015390948f40feacull)
	return false;		// divisible by 193
if ((uint64_t) (n * 0xc96bdb9d3d137e0dull) <= 0x014cab88725af6e7ull)
	return false;		// divisible by 197
if ((uint64_t) (n * 0x2697cc8aef46c0f7ull) <= 0x0149539e3b2d066eull)
	return false;		// divisible by 199
#endif

bool divisibility_sieve(uint64_t t_lo, uint64_t t_hi)
{
	// skip small stooopid requests about 0, 1 and 2 being prime, or composite, or neither, or both.
	if (t_lo > 152 || t_hi != 0) {
		uint64_t v_lo, v_hi;
		// Use the factors of mersenne number 2^60-1 because the modular reduction is simple
		// The factors of 2^60 - 1 are 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321

		// first, reduce mod (2^60-1)
		v_lo = t_lo;
		v_hi = t_hi;
		uint64_t mask60 = (1ull << 60) - 1;
		uint64_t mod_mersenne_60 = v_lo & mask60;
		my_shrd64(&v_hi, &v_lo, 60);
		mod_mersenne_60 += v_lo & mask60;
		my_shrd64(&v_hi, &v_lo, 60);
		mod_mersenne_60 += v_lo;
		// mod_mersenne_60 is a 62 bit number with the same factors than the input number v.
		// i.e. gcd(2^60-1, mod_mersenne_60) == gcd(2^60-1, v)
#if 0
		// these potential factors 3,5,7,11 are already tested tru the wheel sieve.
		if ((uint64_t) (n * 0xaaaaaaaaaaaaaaabull) <= 0x5555555555555555ull)
			return false;	// divisible by 3
		if ((uint64_t) (n * 0xcccccccccccccccdull) <= 0x3333333333333333ull)
			return false;	// divisible by 5
		if ((uint64_t) (n * 0x6db6db6db6db6db7ull) <= 0x2492492492492492ull)
			return false;	// divisible by 7
		if ((uint64_t) (n * 0x2e8ba2e8ba2e8ba3ull) <= 0x1745d1745d1745d1ull)
			return false;	// divisible by 11
#endif
		if ((uint64_t) (mod_mersenne_60 * 0x4ec4ec4ec4ec4ec5ull) <= 0x13b13b13b13b13b1ull)
			return false;	// divisible by 13
		if ((uint64_t) (mod_mersenne_60 * 0xef7bdef7bdef7bdfull) <= 0x0842108421084210ull)
			return false;	// divisible by 31
		if ((uint64_t) (mod_mersenne_60 * 0x8f9c18f9c18f9c19ull) <= 0x063e7063e7063e70ull)
			return false;	// divisible by 41
		if ((uint64_t) (mod_mersenne_60 * 0x4fbcda3ac10c9715ull) <= 0x04325c53ef368eb0ull)
			return false;	// divisible by 61
		if ((uint64_t) (mod_mersenne_60 * 0x6fe4dfc9bf937f27ull) <= 0x01b2036406c80d90ull)
			return false;	// divisible by 151

		// first, reduce mod (2^56-1)
		// 2^56 - 1 = 3 * 5 * 17 * 29 * 43 * 113 * 127 * 15790321
		v_lo = t_lo;
		v_hi = t_hi;
		uint64_t mask56 = (1ull << 56) - 1;
		uint64_t mod_mersenne_56 = v_lo & mask56;
		my_shrd64(&v_hi, &v_lo, 56);
		mod_mersenne_56 += v_lo & mask56;
		my_shrd64(&v_hi, &v_lo, 56);
		mod_mersenne_56 += v_lo;
		// mod_mersenne_56 is a 58 bit number with the same factors than the input number v.
		if ((uint64_t) (mod_mersenne_56 * 0xf0f0f0f0f0f0f0f1ull) <= 0x0f0f0f0f0f0f0f0full)
			return false;	// divisible by 17
		if ((uint64_t) (mod_mersenne_56 * 0x34f72c234f72c235ull) <= 0x08d3dcb08d3dcb08ull)
			return false;	// divisible by 29
		if ((uint64_t) (mod_mersenne_56 * 0x82fa0be82fa0be83ull) <= 0x05f417d05f417d05ull)
			return false;	// divisible by 43
		if ((uint64_t) (mod_mersenne_56 * 0x90fdbc090fdbc091ull) <= 0x0243f6f0243f6f02ull)
			return false;	// divisible by 113
		if ((uint64_t) (mod_mersenne_56 * 0x7efdfbf7efdfbf7full) <= 0x0204081020408102ull)
			return false;	// divisible by 127

		// first, reduce mod (2^36-1)
		// 2^36-1 = 3^3 * 5 * 7 * 13 * 19 * 37 * 73 * 109
		v_lo = t_lo;
		v_hi = t_hi;
		uint64_t mask36 = (1ull << 36) - 1;
		uint64_t mod_mersenne_36 = v_lo & mask36;
		my_shrd64(&v_hi, &v_lo, 36);
		mod_mersenne_36 += v_lo & mask36;
		my_shrd64(&v_hi, &v_lo, 36);
		mod_mersenne_36 += v_lo;
		// mod_mersenne_36 is a 64 bit number where constant modulus are quite cheap 
		if ((uint64_t) (mod_mersenne_36 * 0x86bca1af286bca1bull) <= 0x0d79435e50d79435ull)
			return false;	// divisible by 19
		if ((uint64_t) (mod_mersenne_36 * 0x14c1bacf914c1badull) <= 0x06eb3e45306eb3e4ull)
			return false;	// divisible by 37
		if ((uint64_t) (mod_mersenne_36 * 0x7e3f1f8fc7e3f1f9ull) <= 0x0381c0e070381c0eull)
			return false;	// divisible by 73
		if ((uint64_t) (mod_mersenne_36 * 0xa6c0964fda6c0965ull) <= 0x02593f69b02593f6ull)
			return false;	// divisible by 109

		// first, reduce mod (2^44-1)
		// 2^44-1 : 3 * 5 * 23 * 89 * 397 * 683 * 2113
		v_lo = t_lo;
		v_hi = t_hi;
		uint64_t mask44 = (1ull << 44) - 1;
		uint64_t mod_mersenne_44 = v_lo & mask44;
		my_shrd64(&v_hi, &v_lo, 44);
		mod_mersenne_44 += v_lo & mask44;
		my_shrd64(&v_hi, &v_lo, 44);
		mod_mersenne_44 += v_lo;
		// mod_mersenne_44 is a 64 bit number where constant modulus are quite cheap 
		if ((uint64_t) (mod_mersenne_44 * 0xd37a6f4de9bd37a7ull) <= 0x0b21642c8590b216ull)
			return false;	// divisible by 23
		if ((uint64_t) (mod_mersenne_44 * 0xf47e8fd1fa3f47e9ull) <= 0x02e05c0b81702e05ull)
			return false;	// divisible by 89

		return true;	// might be prime
	} else {
		// t_hi == 0 && t_lo <= 152
		bool stooopid_table[] = {
		true, true, true, true, false, true, false, true, false, false, false, true, false, true, false, false,
		    false, true, false, true, false, false, false, true, false, false, false, false, false, true, false,
		    true, false, false, false, false, false, true, false, false, false, true, false, true, false, false,
		    false, true, false, false, false, false, false, true, false, false, false, false, false, true, false,
		    true, false, false, false, false, false, true, false, false, false, true, false, true, false, false,
		    false, false, false, true, false, false, false, true, false, false, false, false, false, true, false,
		    false, false, false, false, false, false, true, false, false, false, true, false, true, false, false,
		    false, true, false, true, false, false, false, true, false, false, false, false, false, false, false,
		    false, false, false, false, false, false, true, false, false, false, true, false, false, false, false,
		    false, true, false, true, false, false, false, false, false, false, false, false, false, true, false,
		    true, false, false, false, false, false, true, false, false, false, false, false, true, false, false,
		    false, true, false, false, false, false, false, true, false, false, false, false, false, true, false,
		    true, false, false, false, false, false, false, false, false, false, true, false, true, false, false,
		    false, true, false, true, false
		};
		return stooopid_table[t_lo];
	}
}
