// -------------------------------------------------------------------------
//
// Fermat, Euler and Sprp tests
//
// 128 bit utilities
//
// -------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "m128_utils.h"
#include "m_reg.h"

// convert a string to 128 bit integer
//
// string format
//
// 0xa1b2c3d4 .... : hexadecimal number 2712847316
// 1234567890 .... : decimal number 1234567890
// 1e12 .......... : decimal number 1000000000000
// 1.5e12 ........ : decimal number 1500000000000
// 0b010001 ...... : binary number 17
uint128_t convert128(const char *buffer)
{
	char c;
	const char *str = buffer;
	bool dot = false;
	int dotCount = 0;
	unsigned radix = 10;
	uint128_t res = 0;
	// skip leading spaces, if any
	while (*str && isspace(*str)) {
		str++;
	}
	// number prefix
	if (!memcmp(str, "0x", 2) || !memcmp(str, "0X", 2)) {
		radix = 16;
		str += 2;
	}
	if (!memcmp(str, "0b", 2) || !memcmp(str, "0B", 2)) {
		radix = 2;
		str += 2;
	}
	// iterate over the number
	while ((c = *str++) != 0) {
		if (radix == 10 && (c == 'E' || c == 'e')) {
			char *s;
			long exponent = strtol(str, &s, 10);
			while (exponent-- > 0)
				res *= radix;
			str = s + 1;
			break;
		}
		if (c == '.') {
			dot = true;
			continue;
		}
		dotCount += dot ? 1 : 0;
		res *= radix;
		// process digits
		if (c >= '0' && c <= '1') {
			res += c - '0';
		} else if (c >= '2' && c <= '9') {
			if (radix < 10) {
				break;
			}
			res += c - '0';
		} else if (c >= 'A' && c <= 'F') {
			if (radix < 16) {
				break;
			}
			res += c - 'A' + 10;
		} else if (c >= 'a' && c <= 'f') {
			if (radix < 16) {
				break;
			}
			res += c - 'a' + 10;
		} else {
			// unknown digit
			break;
		}
	}
	while (dotCount--) {
		res /= radix;
	}
	// end of input string
	str -= 1;
	if (!str || (!*str) || isspace(*str)) {
		return res;
	} else {
		// not at the end of input string ? bail out. Enough is enough.
		unsigned c = *str;
		printf("Unable to parse the number %s at char %c (0x%2.2x) position %lu\n", buffer, (char)c, c,
		       (unsigned long)(str - buffer));
		fflush(stdout);
		assert(c != 0);
		abort();
	}

	return 0;
}

char *display128(char *temp, uint128_t v, int radix, int digits)
{
	char *pt = temp + 132;
	*--pt = 0;
	char *start_pt = pt;
	if (v == 0) {
		*--pt = '0';
	} else {
		while (v) {
			char c = v % radix;
			v /= radix;
			if (c >= 10) {
				c = c - 10 + 'a';
			} else {
				c = c + '0';
			}
			*--pt = c;
		}
	}
	while (start_pt - pt < digits) {
		*--pt = '0';
	}

	if (radix == 2) {
		*--pt = 'b';
		*--pt = '0';
	} else if (radix == 16) {
		*--pt = 'x';
		*--pt = '0';
	}
	return pt;
}

void my_printf(uint128_t v)
{
	char temp[140];
	printf("%s", display128(temp, v, 16, 32));
}

void my_printf(uint64_t v_lo, uint64_t v_hi)
{
	char temp[140];
	uint128_t v = ((uint128_t) v_hi << 64) + v_lo;
	printf("%s", display128(temp, v, 16, 32));
}

uint64_t my_clz128(uint128_t v)
{
	return my_clz128((uint64_t) v, (uint64_t) (v >> 64));
}

uint64_t my_ctz128(uint128_t v)
{
	return my_ctz128((uint64_t) v, (uint64_t) (v >> 64));
}
