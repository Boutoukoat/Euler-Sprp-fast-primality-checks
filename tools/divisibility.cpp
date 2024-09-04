
#include <stdint.h>
#include <m_reg.h>

// ------------------------------------------------------------------------
// fast divisibility of 128 bits numbers by small constants 13, 31, 41, 61 .... 
//
// // assuming that divisibility by 2,3,5,11 has already been done
//
// This is not the fastest way, modular operations on constants modulus have many tricks (TODO)
// ------------------------------------------------------------------------
bool divisibility_sieve(uint64_t t_lo, uint64_t t_hi)
{
	// skip small stooopid requests about 0, 1 and 2 being prime, or composite, or neither, or both.
	if (t_lo > 128 || t_hi != 0)
	{
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
	// mod_mersenne_60 is now a 62 bit number with the same factors than the input number v.
	// i.e. gcd(2^60-1, mod_mersenne_60) == gcd(2^60-1, v)

	if (mod_mersenne_60 % 13 == 0)
		return false;	// composite for sure
	if (mod_mersenne_60 % 31 == 0)
		return false;	// composite for sure
	if (mod_mersenne_60 % 41 == 0)
		return false;	// composite for sure
	if (mod_mersenne_60 % 61 == 0)
		return false;	// composite for sure

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
	// mod_mersenne_56 is now a 58 bit number with the same factors than the input number v.
	// mod_mersenne_36 is now a 64 bit number, where constant modulus are quite cheap 
	if (mod_mersenne_56 % 17 == 0)
		return false;	// composite for sure
	if (mod_mersenne_56 % 29 == 0)
		return false;	// composite for sure
	if (mod_mersenne_56 % 43 == 0)
		return false;	// composite for sure
	if (mod_mersenne_56 % 113 == 0)
		return false;	// composite for sure
	if (mod_mersenne_56 % 127 == 0)
		return false;	// composite for sure

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
	// mod_mersenne_36 is now a 64 bit number where constant modulus are quite cheap 
	if (mod_mersenne_36 % 19 == 0)
		return false;	// composite for sure
	if (mod_mersenne_36 % 37 == 0)
		return false;	// composite for sure
	if (mod_mersenne_36 % 73 == 0)
		return false;	// composite for sure
	if (mod_mersenne_36 % 109 == 0)
		return false;	// composite for sure

	return true;		// might be prime
	}
	else
	{

	bool stooopid_table[] =
	{true, true, true, true, false, true, false, true, false, false, false, true, false, true, false, false, false, true, false, true, false, false, false, true, false, false, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, true, false, true, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, true, false, true, false, false, false, false, false, true, false, false, false, true, false, false, false, false, false, true, false, false, false, false, false, false, false, true, false, false, false, true, false, true, false, false, false, true, false, true, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false};
	return stooopid_table[t_lo];
	}
}


