
#include <stdint.h>
#include <m_reg.h>

// ------------------------------------------------------------------------
// fast divisibility of 128 bits numbers by small constants 13, 31, 41, 61 .... 
//
// This is not the fastest way, modular operations on constants modulus have many tricks (TODO)
// ------------------------------------------------------------------------
bool divisibility_sieve(uint64_t t_lo, uint64_t t_hi)
{
	uint64_t v_lo, v_hi;
	// Use the factors of mersenne number 2^60-1 because the modular reduction is simple
	// The factors of 2^60 - 1 are 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321

	// first, reduce mod (2^60-1)
	v_lo = t_lo;
	v_hi = t_hi;
	uint64_t mask60 = (1ull << 60) - 1;
	uint64_t mod60 = v_lo & mask60;
	my_shrd64(&v_hi, &v_lo, 60);
	mod60 += v_lo & mask60;
	my_shrd64(&v_hi, &v_lo, 60);
	mod60 += v_lo;
	// mod60 is now a 62 bit number with the same factors than the input number v.
	// i.e. gcd(2^60-1, mod60) == gcd(2^60-1, v)

	// since factors 3, 2, 5, 11 are skipped, let's start with 13
	if (mod60 % 13 == 0 && (t_hi > 0 || t_lo > 13))
		return false;	// composite for sure
	if (mod60 % 31 == 0 && (t_hi > 0 || t_lo > 31))
		return false;	// composite for sure
	if (mod60 % 41 == 0 && (t_hi > 0 || t_lo > 41))
		return false;	// composite for sure
	if (mod60 % 61 == 0 && (t_hi > 0 || t_lo > 61))
		return false;	// composite for sure

	// first, reduce mod (2^56-1)
	// 2^56 - 1 = 3 * 5 * 17 * 29 * 43 * 113 * 127 * 15790321
	v_lo = t_lo;
	v_hi = t_hi;
	uint64_t mask56 = (1ull << 56) - 1;
	uint64_t mod56 = v_lo & mask56;
	my_shrd64(&v_hi, &v_lo, 56);
	mod56 += v_lo & mask56;
	my_shrd64(&v_hi, &v_lo, 56);
	mod56 += v_lo;
	// mod56 is now a 58 bit number with the same factors than the input number v.
	if (mod56 % 17 == 0 && (t_hi > 0 || t_lo > 17))
		return false;	// composite for sure
	if (mod56 % 29 == 0 && (t_hi > 0 || t_lo > 29))
		return false;	// composite for sure
	if (mod56 % 43 == 0 && (t_hi > 0 || t_lo > 43))
		return false;	// composite for sure

	return true;		// might be prime
}
