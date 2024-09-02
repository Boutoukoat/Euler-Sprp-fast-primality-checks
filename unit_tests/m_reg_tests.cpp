
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#include "m_reg.h"

int main(int argc, char **argv)
{
	uint64_t check_count = 0;
	bool verbose = false;

	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-v")) {
			verbose = true;
		}
	}

#if defined(INLINE_ASM)
#if INLINE_ASM==0
	printf("INLINE_ASM is 0 (slow ....)\n");
#else
	// inline assembly language is enabled
#endif
#else
	printf("INLINE_ASM is undefined (slow ....)\n");
#endif

	uint64_t v[] = { 0, 1, 2, (uint64_t) - 1, (uint64_t) 1 << 63, (uint64_t) - 1 << 1, (uint64_t) - 1 >> 1 };

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < sizeof(v) / sizeof(v[0]); j++) {
			for (unsigned c = 0; c < 2; c++) {
				uint128_t t = (uint128_t) v[j] + (uint128_t) v[i] + (uint128_t) c;
				uint64_t r = 0xc0c0c0c0c0c0c0;
				uint8_t rc = 0xc0;
				uint64_t e = (uint64_t) t;
				uint8_t ec = (t >> 64) ? 1 : 0;

				r = 0xc0c0c0c0c0c0c0;
				rc = 0xc0;
				rc = my_adc64(c, v[i], v[j], &r);
				assert(r == e);
				check_count++;
				assert(rc == ec);
				check_count++;
			}
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < sizeof(v) / sizeof(v[0]); j++) {
			for (unsigned c = 0; c < 2; c++) {
				uint128_t t = (uint128_t) v[i] - (uint128_t) v[j] - (uint128_t) c;
				uint64_t r = 0xc0c0c0c0c0c0c0;
				uint8_t rc = 0xc0;
				uint64_t e = (uint64_t) t;
				uint8_t ec = (t >> 64) ? 1 : 0;

				r = 0xc0c0c0c0c0c0c0;
				rc = 0xc0;
				rc = my_sbb64(c, v[i], v[j], &r);
				assert(r == e);
				check_count++;
				assert(rc == ec);
				check_count++;
			}
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < sizeof(v) / sizeof(v[0]); j++) {
			uint128_t t = (uint128_t) v[i] * (uint128_t) v[j];
			uint64_t e_lo = (uint64_t) t;
			uint64_t e_hi = (uint64_t) (t >> 64);

			uint64_t r_lo = 0xc0c0c0c0c0c0c0;
			uint64_t r_hi = 0xc0c0c0c0c0c0c0;
			r_lo = my_mul64(v[i], v[j], &r_hi);
			assert(r_hi == e_hi);
			check_count++;
			assert(r_lo == e_lo);
			check_count++;
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < 4; j++) {
			uint128_t t = (uint128_t) v[i] << j;
			uint64_t e_lo = (uint64_t) t;
			uint64_t e_hi = (uint64_t) (t >> 64);

			uint64_t r_lo = 0xc0c0c0c0c0c0c0;
			r_lo = my_shl64(v[i], j);
			assert(r_lo == e_lo);
			check_count++;
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < 4; j++) {
			uint128_t t = (uint128_t) v[i] >> j;
			uint64_t e_lo = (uint64_t) t;
			uint64_t e_hi = (uint64_t) (t >> 64);

			uint64_t r_lo = 0xc0c0c0c0c0c0c0;
			r_lo = my_shr64(v[i], j);
			assert(r_lo == e_lo);
			check_count++;
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < sizeof(v) / sizeof(v[0]); j++) {
			for (unsigned k = 0; k < 4; k++) {
				uint128_t t = (((uint128_t) v[i] << 64) + v[j]) << k;
				uint64_t e_lo = (uint64_t) t;
				uint64_t e_hi = (uint64_t) (t >> 64);
				uint64_t t1 = v[i];
				uint64_t t0 = v[j];

				my_shld64(&t1, &t0, k);
				assert(t0 == e_lo);
				check_count++;
				assert(t1 == e_hi);
				check_count++;
			}
		}
	}

	for (unsigned i = 0; i < sizeof(v) / sizeof(v[0]); i++) {
		for (unsigned j = 0; j < sizeof(v) / sizeof(v[0]); j++) {
			for (unsigned k = 0; k < 4; k++) {
				uint128_t t = (((uint128_t) v[i] << 64) + v[j]) >> k;
				uint64_t e_lo = (uint64_t) t;
				uint64_t e_hi = (uint64_t) (t >> 64);
				uint64_t t1 = v[i];
				uint64_t t0 = v[j];

				my_shrd64(&t1, &t0, k);
				assert(t0 == e_lo);
				check_count++;
				assert(t1 == e_hi);
				check_count++;
			}
		}
	}

	uint64_t x2 = (uint64_t)1 << 50;
	assert(my_clz64(x2) == 13);
				check_count++;
	assert(my_ctz64(x2) == 50);
				check_count++;
	uint32_t x1 = (uint32_t)1 << 20;
	assert(my_clz32(x1) == 11);
				check_count++;
	assert(my_ctz32(x1) == 20);
				check_count++;

	assert(check_count != 0);
	printf("%lu checks passed\n", check_count);
	return 0;
}
