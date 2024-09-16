#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "m128_utils.h"
#include "avx2_sprp.cpp"

// min and max primes for every bit range

struct {
	const char *p_min;
	const char *p_max;
} primes_min_max[] = {
	{"0x3", "0x3"},
	{"0x5", "0x7"},
	{"0xB", "0xD"},
	{"0x11", "0x1F"},
	{"0x25", "0x3D"},
	{"0x43", "0x7F"},
	{"0x83", "0xFB"},
	{"0x101", "0x1FD"},
	{"0x209", "0x3FD"},
	{"0x407", "0x7F7"},
	{"0x805", "0xFFD"},
	{"0x1003", "0x1FFF"},
	{"0x2011", "0x3FFD"},
	{"0x401B", "0x7FED"},
	{"0x8003", "0xFFF1"},
	{"0x10001", "0x1FFFF"},
	{"0x2001D", "0x3FFFB"},
	{"0x40003", "0x7FFFF"},
	{"0x80015", "0xFFFFD"},
	{"0x100007", "0x1FFFF7"},
	{"0x200011", "0x3FFFFD"},
	{"0x40000F", "0x7FFFF1"},
	{"0x800009", "0xFFFFFD"},
	{"0x100002B", "0x1FFFFD9"},
	{"0x2000023", "0x3FFFFFB"},
	{"0x400000F", "0x7FFFFD9"},
	{"0x800001D", "0xFFFFFC7"},
	{"0x10000003", "0x1FFFFFFD"},
	{"0x2000000B", "0x3FFFFFDD"},
	{"0x40000003", "0x7FFFFFFF"},
	{"0x8000000B", "0xFFFFFFFB"},
	{"0x10000000F", "0x1FFFFFFF7"},
	{"0x200000011", "0x3FFFFFFD7"},
	{"0x400000019", "0x7FFFFFFE1"},
	{"0x800000035", "0xFFFFFFFFB"},
	{"0x100000001F", "0x1FFFFFFFE7"},
	{"0x2000000009", "0x3FFFFFFFD3"},
	{"0x4000000007", "0x7FFFFFFFF9"},
	{"0x8000000017", "0xFFFFFFFFA9"},
	{"0x1000000000F", "0x1FFFFFFFFEB"},
	{"0x2000000001B", "0x3FFFFFFFFF5"},
	{"0x4000000000F", "0x7FFFFFFFFC7"},
	{"0x8000000001D", "0xFFFFFFFFFEF"},
	{"0x100000000007", "0x1FFFFFFFFFC9"},
	{"0x20000000003B", "0x3FFFFFFFFFEB"},
	{"0x40000000000F", "0x7FFFFFFFFF8D"},
	{"0x800000000005", "0xFFFFFFFFFFC5"},
	{"0x1000000000015", "0x1FFFFFFFFFFAF"},
	{"0x2000000000045", "0x3FFFFFFFFFFE5"},
	{"0x4000000000037", "0x7FFFFFFFFFF7F"},
	{"0x8000000000015", "0xFFFFFFFFFFFD1"},
	{"0x10000000000015", "0x1FFFFFFFFFFF91"},
	{"0x20000000000005", "0x3FFFFFFFFFFFDF"},
	{"0x4000000000009F", "0x7FFFFFFFFFFFC9"},
	{"0x80000000000003", "0xFFFFFFFFFFFFFB"},
	{"0x100000000000051", "0x1FFFFFFFFFFFFF3"},
	{"0x200000000000009", "0x3FFFFFFFFFFFFE5"},
	{"0x400000000000045", "0x7FFFFFFFFFFFFC9"},
	{"0x800000000000083", "0xFFFFFFFFFFFFFA3"},
	{"0x1000000000000021", "0x1FFFFFFFFFFFFFFF"},
	{"0x200000000000000F", "0x3FFFFFFFFFFFFFC7"},
	{"0x4000000000000087", "0x7FFFFFFFFFFFFFE7"},
	{"0x800000000000001D", "0xFFFFFFFFFFFFFFC5"},
	{"0x1000000000000000D", "0x1FFFFFFFFFFFFFFCF"},
	{"0x20000000000000083", "0x3FFFFFFFFFFFFFFFB"},
	{"0x40000000000000009", "0x7FFFFFFFFFFFFFFED"},
	{"0x80000000000000003", "0xFFFFFFFFFFFFFFFE9"},
	{"0x100000000000000021", "0x1FFFFFFFFFFFFFFFED"},
	{"0x20000000000000001D", "0x3FFFFFFFFFFFFFFFDD"},
	{"0x400000000000000019", "0x7FFFFFFFFFFFFFFF19"},
	{"0x80000000000000000B", "0xFFFFFFFFFFFFFFFFA3"},
	{"0x100000000000000000F", "0x1FFFFFFFFFFFFFFFFBB"},
	{"0x200000000000000001D", "0x3FFFFFFFFFFFFFFFFDD"},
	{"0x4000000000000000025", "0x7FFFFFFFFFFFFFFFF9F"},
	{"0x8000000000000000021", "0xFFFFFFFFFFFFFFFFFF1"},
	{"0x1000000000000000000F", "0x1FFFFFFFFFFFFFFFFFDF"},
	{"0x2000000000000000000B", "0x3FFFFFFFFFFFFFFFFFF5"},
	{"0x40000000000000000007", "0x7FFFFFFFFFFFFFFFFFBD"},
	{"0x80000000000000000017", "0xFFFFFFFFFFFFFFFFFFBF"},
	{"0x10000000000000000000D", "0x1FFFFFFFFFFFFFFFFFFCD"},
	{"0x200000000000000000011", "0x3FFFFFFFFFFFFFFFFFFC7"},
	{"0x400000000000000000009", "0x7FFFFFFFFFFFFFFFFFFC9"},
	{"0x80000000000000000004B", "0xFFFFFFFFFFFFFFFFFFFDD"},
	{"0x1000000000000000000003", "0x1FFFFFFFFFFFFFFFFFFFED"},
	{"0x20000000000000000000AB", "0x3FFFFFFFFFFFFFFFFFFFDD"},
	{"0x400000000000000000001B", "0x7FFFFFFFFFFFFFFFFFFFBD"},
	{"0x8000000000000000000027", "0xFFFFFFFFFFFFFFFFFFFED5"},
	{"0x10000000000000000000007", "0x1FFFFFFFFFFFFFFFFFFFFFF"},
	{"0x2000000000000000000001D", "0x3FFFFFFFFFFFFFFFFFFFFDF"},
	{"0x40000000000000000000085", "0x7FFFFFFFFFFFFFFFFFFFFD3"},
	{"0x8000000000000000000003B", "0xFFFFFFFFFFFFFFFFFFFFFAD"},
	{"0x100000000000000000000019", "0x1FFFFFFFFFFFFFFFFFFFFFE7"},
	{"0x200000000000000000000069", "0x3FFFFFFFFFFFFFFFFFFFFFFD"},
	{"0x400000000000000000000081", "0x7FFFFFFFFFFFFFFFFFFFFFF1"},
	{"0x800000000000000000000009", "0xFFFFFFFFFFFFFFFFFFFFFFEF"},
	{"0x100000000000000000000003D", "0x1FFFFFFFFFFFFFFFFFFFFFF73"},
	{"0x2000000000000000000000069", "0x3FFFFFFFFFFFFFFFFFFFFFFCD"},
	{"0x4000000000000000000000007", "0x7FFFFFFFFFFFFFFFFFFFFFF8D"},
	{"0x80000000000000000000000FF", "0xFFFFFFFFFFFFFFFFFFFFFFFF1"},
	{"0x10000000000000000000000115", "0x1FFFFFFFFFFFFFFFFFFFFFFFBB"},
	{"0x20000000000000000000000051", "0x3FFFFFFFFFFFFFFFFFFFFFFFDF"},
	{"0x4000000000000000000000010B", "0x7FFFFFFFFFFFFFFFFFFFFFFF9F"},
	{"0x80000000000000000000000051", "0xFFFFFFFFFFFFFFFFFFFFFFFFEF"},
	{"0x10000000000000000000000006F", "0x1FFFFFFFFFFFFFFFFFFFFFFFFF3"},
	{"0x200000000000000000000000027", "0x3FFFFFFFFFFFFFFFFFFFFFFFF8B"},
	{"0x400000000000000000000000063", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFF"},
	{"0x800000000000000000000000027", "0xFFFFFFFFFFFFFFFFFFFFFFFFFC5"},
	{"0x1000000000000000000000000021", "0x1FFFFFFFFFFFFFFFFFFFFFFFFFE1"},
	{"0x2000000000000000000000000093", "0x3FFFFFFFFFFFFFFFFFFFFFFFFFEB"},
	{"0x400000000000000000000000001B", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFDB"},
	{"0x8000000000000000000000000033", "0xFFFFFFFFFFFFFFFFFFFFFFFFFFB5"},
	{"0x10000000000000000000000000019", "0x1FFFFFFFFFFFFFFFFFFFFFFFFFF7B"},
	{"0x20000000000000000000000000119", "0x3FFFFFFFFFFFFFFFFFFFFFFFFFFF5"},
	{"0x4000000000000000000000000002B", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFFBD"},
	{"0x80000000000000000000000000047", "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFD"},
	{"0x100000000000000000000000000021", "0x1FFFFFFFFFFFFFFFFFFFFFFFFFFEE9"},
	{"0x20000000000000000000000000001D", "0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFB"},
	{"0x400000000000000000000000000019", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFBB"},
	{"0x800000000000000000000000000009", "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF89"},
	{"0x10000000000000000000000000001C3", "0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFB7"},
	{"0x2000000000000000000000000000029", "0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFFD"},
	{"0x4000000000000000000000000000115", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFBD"},
	{"0x80000000000000000000000000000A5", "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC5"},
	{"0x10000000000000000000000000000043", "0x1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF7"},
	{"0x2000000000000000000000000000001B", "0x3FFFFFFFFFFFFFFFFFFFFFFFFFFFFF77"},
	{"0x40000000000000000000000000000007", "0x7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"},
	{"0x8000000000000000000000000000001D", "0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF61"}
};

const char *pseudoprimes[] = {
	"0x4774c39",
	"0x26d55381",
	"0x108de3061",
	"0x2d7843ff9",
	"0xab40d0929",
	"0x6c18494e99",
	"0xab1363fb79",
	"0x58bd9ad4f69",
	"0x18cd7c8c3159",
	"0x46b704509d11",
	"0x19c8847f327e1",
	"0x2e69f8daf2b79",
	"0xd245d34ea1199",
	"0x3609e7f2b03399",
	"0x8e8e5f5fca61e9",
	"0x46750bf64ef68b9",
	"0x12980948d17c7a39",
	"0x720bb118c7fc4471",
	"0x18e9c4fdbdbd65d99",
	"0x2494a86f64d8c1be9",
	"0x14f4cbe19cee228c01",
	"0x35570c9f5febbee731",
	"0x1c9c2d7c5383a810339",
	"0x2b22b5b362c483faa01",
	"0x83c0f1fdd208f20f059",
	"0x2a9f0bf59e36a31478a9",
	"0xc4f69ea753c651a0bfc9",
	"0x21bb0e110f3eabcd2d179",
	"0xb892633c0c277254033d9",
	"0x4a0dd0ae918f26ef758cb1",
	"0x83e8abbabe4a03cccd9ed9",
	"0x25f2e32a94c0474829398f1",
	"0x9749ca2234a3afb9358d2b1",
	"0x6aa8ade3e1102293b114fde1",
	"0xec01fee44a0c199a651d7f89",
	"0x413eb180e9aabf94bc8acc949",
	"0x1dd2612de0550e35d6514209d1",
	"0x314cb7b1ff163f039094a58ed1",
	"0xfa4bc7a0b0b8b0e8cc309ba9c1",
	"0x35de61aa9d961a0bcbd66ea95a9",
	"0x1c476bd48b340529d3e46e28c4f1",
	"0x4e69379f154a76c1967e2b07a6a1",
	"0x1e1d5c87aa298c26e3f12fc4b7d81",
	"0x45f4143df35f3684938bfb3cc1c21",
	"0x98846af16c6f15f0b1380fcdefbd1",
	"0x7df5320417d531ac139b1314cd9049",
	"0x15429dd1389604d6bb99b46d698ec99",
	"0x62fe7df4a56bbffae277136732e7969",
	"0x12472edeeb3ece07631257ba49329d51",
	"0x2a236fe00c37677ec1a9a9dace193089",
	"0xd2a24bb1346b748acb34bf69d9cf6d19",
	"0x53c3e9",
	"0x8524f9",
	"0xd0f729",
	"0x7085a01",
	"0x8ffe7d1",
	"0x14ecbb41",
	"0x3a5fe551",
	"0x62636929",
	"0xd6e23e09",
	"0x1c1a8d439",
	"0x2431350c1",
	"0x759123211",
	"0xe6e564b79",
	"0x1b02a5bed9",
	"0x2abece5e81",
	"0x6603229ef1",
	"0xeb773375a9",
	"0x1ed9d37f9d9",
	"0x3f940b02881",
	"0x4431800da21",
	"0x9065ff8b8d1",
	"0x198f56176a71",
	"0x30e87de2c9f9",
	"0x4471bb7a1ce9",
	"0x89d59706fa81",
	"0x1f6f03e167ce9",
	"0x329ecd0c06449",
	"0x70fee10c07321",
	"0xdd12f431eb7d9",
	"0x17257ab0474549",
	"0x265ad224c38b91",
	"0x680689f0e9bb09",
	"0xf3c40550896ff9",
	"0x135107abac0e9b9",
	"0x242edee2a4b4621",
	"0x5f84ed7f36b5ce1",
	"0x8ddb53df35bdef9",
	"0x1d389f972ecbcb51",
	"0x31fe1b8785a549b9",
	"0x415b8ac9cd51bff1",
	"0xb3123e2f48d0b739",
	"0x11a69d5b80e61e949",
	"0x39352a7737709f781",
	"0x7e599012cb7e9f819",
	"0xea9cd031028a2bb39",
	"0x17326ebc39a16b6cd9",
	"0x3933bb9f4abaeb6e31",
	"0x41e277624c473b5539",
	"0xc4182bca32c561a5d1",
	"0x112833e9740e3e68e61",
	"0x2e4e80142a721427e61",
	"0x567673d9e88d0572ec9",
	"0xc0e9ba151dd5186d989",
	"0x172bbc76705cd6c21729",
	"0x2db0340db246528f65a1",
	"0x660cd3242e45b5ee9129",
	"0xf6b590e39094da160b61",
	"0x10999f885c91a71e91689",
	"0x2ddf3a32ea03fcb747289",
	"0x5068dc219afafaef00c81",
	"0xc5f5abb538da4ad61c5d9",
	"0x1b5e3dc979dd43ed0244b1",
	"0x36a0aaa15078ca1d4df6d9",
	"0x538ac459f48af0dc680cd9",
	"0xe8761e8aed6ba2317410e9",
	"0x10276fd326b7d43fd9299f1",
	"0x25f3b9b803eb19957911629",
	"0x54c936c14bc176a593c9889",
	"0xd1055dba5cda294d38bba09",
	"0x1ae9ccad8b5452a6edccff91",
	"0x22df0a1751349def6bcbb619",
	"0x62313861c3af8d2fe5c0bb11",
	"0x96c9f73575b61d11cd8f8119",
	"0x1eb329ad43e31815db1a63749",
	"0x2b02c5d5bae650d7b60893271",
	"0x72a55a358ca9a09089e5099d9",
	"0x8f66590d84740150d93855f89",
	"0x12564a342e9e8b130664b212d9",
	"0x3ee440e8a9a49c3246983af0c9",
	"0x7b1108589767baa1f84e74b1e1",
	"0xab97550cd47f02e283a46059b9",
	"0x17e3541a042e297678944988659",
	"0x273558952da5d0ee07aed883201",
	"0x5e780622e5a4add5240b7b04141",
	"0xbc782354da415be88b8f35cd861",
	"0x16088dce02ba122050dfbce535e9",
	"0x28d1fce8b2f608f81c96866144f1",
	"0x54e5fa2fbbfe438be496be611e29",
	"0x8045b489db349ef10cbaeedb3649",
	"0x10e729b5bbf6e5663c828af8c40c9",
	"0x2d3836e2b27e6df2acf7c7ec2a031",
	"0x5b0f7733a3f31ad501bb6b93b7b31",
	"0xdece847776589255c239276626d59",
	"0x106d2362df56dad04df37ba49e2d49",
	"0x3210b5ac2ce54a90065b0a20f23bf1",
	"0x503eb4438d5989125bbfab80886961",
	"0x96403e0f207f72be50ac2344af7ce1",
	"0x1622ae9bf4996323a464e926cb2dc39",
	"0x2d8d25d98357f8c01758df57fb35a91",
	"0x490f93ecafcd398a1f19f8b5568d181",
	"0xd7c650e77c739b80bb1ac43bc070959",
	"0x1d62c36f064b898b1fe069aa07974be9",
	"0x27aef682f2317bd1548bb3a40879ab71",
	"0x6d4b3fd2cc046fbb2762cc0e31981f29",
	"0x87b9f48eb0ef03306f4273ae919850e9"
};

static const char *square_root[] = {
	"0x3201",
	"0x55001",
	"0x80901",
	"0xd0201",
	"0x21c501",
	"0x37cb01",
	"0x677101",
	"0x1010101",
	"0x1d50e01",
	"0x2325301",
	"0x58b8701",
	"0xdbf2401",
	"0x15ee5001",
	"0x40014001",
	"0x5862e701",
	"0xbdde7201",
	"0x15a8fab01",
	"0x229cce601",
	"0x535f9ca01",
	"0xa63419001",
	"0x1a89771201",
	"0x2142355501",
	"0x655ec35801",
	"0xf43b5cd901",
	"0x1c5366edf01",
	"0x2188cfa8b01"
};

static const char *composites[] = {
	"0x2e5",
	"0x7a1",
	"0x3111",
	"0x977d",
	"0x13629",
	"0x15ca5",
	"0x2e185",
	"0x79b45",
	"0x1023f9",
	"0x14159d",
	"0x3bc019",
	"0x5e7361",
	"0xa30569",
	"0x1eef4ed",
	"0x290d1a1",
	"0x760d7ed",
	"0xd7e6551",
	"0x11c69e05",
	"0x2e343dc9",
	"0x7afb5819",
	"0x91865ead",
	"0x16b7d2605",
	"0x302ba1639",
	"0x7335c6ec1",
	"0xfc30a1b1d",
	"0x191bca20d1",
	"0x2e5674ab6d",
	"0x57c99421d1",
	"0xbbd638af31",
	"0x10212b3ac15",
	"0x2827eeb48f5",
	"0x5a83d5124a5",
	"0xad19d6bace9",
	"0x1639887da3cd",
	"0x225534723d05",
	"0x55e70600ec69",
	"0x955d9869a211",
	"0x19151312a37ad",
	"0x296e29f5e2201",
	"0x7a32b6df1072d",
	"0xc2fcf9e961fb5",
	"0x1252c06aa5bf95",
	"0x27d375e1d4f291",
	"0x417bb10805c459",
	"0xa2294765d46019",
	"0x1e5bfc27ce31a0d",
	"0x2f5575906b26785",
	"0x6bc52943d43a941",
	"0x9cd9a3cf45788a1",
	"0x1fa857a9b8dba011",
	"0x2a4d9e7532d56dd1",
	"0x72e283c5ed0a64ad",
	"0xf10df8dcd7469c0d",
	"0x18ede8f7af4ffed95",
	"0x2c68ab2f0448a0955",
	"0x66cb1007dc005a1bd",
	"0xf06d4a2ed89ddad61",
	"0x1fef3011961223aaed",
	"0x23dadcaf165c013acd",
	"0x720f9cff1bfcf3aca1",
	"0xc4f98c38d98bc42f65",
	"0x196e4d2e12229335099",
	"0x22b24ae7760e007f51d",
	"0x655fce197de60c5f635",
	"0xb64acc572af7a126e61",
	"0x123ad437221c23c412b1",
	"0x2d830926316e36b2b245",
	"0x77d62b25f4810c836c01",
	"0x8db3ec4120ff62349715",
	"0x10661691b3862fa9dfaf9",
	"0x32c8724f5ab1c19d0c9c5",
	"0x7ff4a97d3ab56f43b3875",
	"0x859e80c0879d8a51d9ed5",
	"0x10cf98a3a690b55a0447bd",
	"0x3f9c252ea9dfb7fc787edd",
	"0x56fc1c522659fc74c536ed",
	"0xf7b99c1c6b11219d78f6bd",
	"0x100220a74f94954e0897931",
	"0x3fd252b7bda698b457f7295",
	"0x5f8f080686db659eb22f291",
	"0xae9ed763e3d6ab40ebea921",
	"0x1479fa60f8d55ed446d1a7fd",
	"0x2473f26672829a013ae34ac1",
	"0x441cd23b4748ac7ee98bea49",
	"0xab48d75d961f88e4e5101905",
	"0x1341f100e0a85698ac7774bc1",
	"0x36e4bbd66209388eddaa3ebbd",
	"0x7c9dd349181255f9ec02cb4a5",
	"0xf1028134915b2a19186a78fa9",
	"0x17fab32de5203e48305edfa4cd",
	"0x3c2aecadfd1f786615246e22d5",
	"0x5cd3a34fe89679f4ebd8695b65",
	"0xc283f533dbd01e528dbd7ccf09",
	"0x16e4213b80fa67f5426ded2aa85",
	"0x3c8e354ad86fbb6648addc286d1",
	"0x442fbec498ec80ca157edeb51e1",
	"0xcea3eb0ad91a24dcffc0b8bf819",
	"0x1ea41b2325175fdd63404e837171",
	"0x34b96fa67cb43c8cab084f8b34c9",
	"0x4c21cc5a8d5e279a97c7a4f0ce2d",
	"0x81b50e2499695cf42ced25282165",
	"0x18fa147062ba6bd0d51b87727a799",
	"0x33089a1992a75aeab201a67962d09",
	"0x531bd1ac74b01b2a0d1cbb2c04d95",
	"0xbed2269ffd53812b768c8fc6b9211",
	"0x1e296c9e5c1140e7efad8f88a047f5",
	"0x24e39831fd5e8c69423bce95eac80d",
	"0x428c3611f02e057def2857d96594ad",
	"0xf03340ca1104147db81553c8f73fbd",
	"0x12d4e5b2961a7122fd0f36bae542109",
	"0x299f80e7cf30d1c1358d6a64b773309",
	"0x6f8849f5143b1221800f041429b8fe9",
	"0xdfe24e79d96c09249544735544a7e95",
	"0x1a0f21dca9992a27fe2e14c96adbf311",
	"0x37afde07c5f0e913d5316cd309f5d7d9",
	"0x4018cae1ee1f03e8d082d3c2784f6385",
	"0xd9078d3e1214b2477115a9793bf90621"
};

static uint64_t self_tests(void)
{
	uint64_t check_count = 0;
	bool b;

#ifdef __AVX512F__

	// check endianness , different in C intrinsics vs. gdb display vs. documentation vs. memory layout.
	__m512i g;
	__m256i e;
	__mmask16 k = 0x5555;

	g = _mm512_set_epi32(1,2,3,4,5,6,7,8, 0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf);
	e = _mm512_extracti32x8_epi32(g, 0);

	assert(_mm256_extract_epi32(e, 0) == 0xf);
	check_count++;
	assert(_mm256_extract_epi32(e, 1) == 0xe);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0xd);
	check_count++;
	assert(_mm256_extract_epi32(e, 3) == 0xc);
	check_count++;

	g = _mm512_maskz_shuffle_epi32(k, g, _MM_PERM_DCBA);
	e = _mm512_extracti32x8_epi32(g, 0);

	assert(_mm256_extract_epi32(e, 0) == 0xf);
	check_count++;
	assert(_mm256_extract_epi32(e, 1) == 0x0);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0xd);
	check_count++;
	assert(_mm256_extract_epi32(e, 3) == 0x0);
	check_count++;

	g = _mm512_set_epi64(0xa, 0xb, 0xc, 0xd, 0x1, 0x2, 0x3, 0x4);
	e = _mm512_extracti32x8_epi32(g, 0);
	assert(_mm256_extract_epi32(e, 0) == 0x4);
	check_count++;
	assert(_mm256_extract_epi32(e, 1) == 0x0);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0x3);
	check_count++;
	assert(_mm256_extract_epi32(e, 3) == 0x0);
	check_count++;

	g = _mm512_set_epi32(1,2,3,4,5,6,7,8, 0x8, 0x9, 0xa, 0xb, 0xc, 0xd, 0xe, 0xf);
	g = _mm512_shuffle_epi32(g, _MM_PERM_DDBB);
	e = _mm512_extracti32x8_epi32(g, 0);
	assert(_mm256_extract_epi32(e, 0) == 0xe);
	check_count++;
	assert(_mm256_extract_epi32(e, 1) == 0xe);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0xc);
	check_count++;
	assert(_mm256_extract_epi32(e, 3) == 0xc);
	check_count++;
	assert(_mm256_extract_epi32(e, 4) == 0xa);
	check_count++;
	assert(_mm256_extract_epi32(e, 5) == 0xa);
	check_count++;
	assert(_mm256_extract_epi32(e, 6) == 0x8);
	check_count++;
	assert(_mm256_extract_epi32(e, 7) == 0x8);
	check_count++;
	e = _mm512_extracti32x8_epi32(g, 1);
	assert(_mm256_extract_epi32(e, 0) == 0x7);
	check_count++;
	assert(_mm256_extract_epi32(e, 1) == 0x7);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0x5);
	check_count++;
	assert(_mm256_extract_epi32(e, 3) == 0x5);
	check_count++;
	assert(_mm256_extract_epi32(e, 4) == 0x3);
	check_count++;
	assert(_mm256_extract_epi32(e, 5) == 0x3);
	check_count++;
	assert(_mm256_extract_epi32(e, 6) == 0x1);
	check_count++;
	assert(_mm256_extract_epi32(e, 7) == 0x1);
	check_count++;

	// check sprp exit conditions
	uint64_t mask;
	__m512i ma[1], mb[1];
	mask = 0;
	ma[0] = _mm512_set_epi32(0,1,0,2,0,3,0,4, 0, 1, 0, 2, 0, 3, 0, 4);
	mb[0] = _mm512_set_epi32(0,8,0,8,0,8,0,8, 0, 8, 0, 8, 0, 8, 0, 8);
	b = avx512_cmp_next1(&mask, ma, mb);
	assert(b == false);
	check_count++;
	assert(mask == 0x00);
	check_count++;
	ma[0] = _mm512_set_epi32(0,1,0,2,0,3,0,4, 0, 1, 0, 8, 0, 3, 0, 4);
	b = avx512_cmp_next1(&mask, ma, mb);
	assert(b == false);
	check_count++;
	assert(mask == 0x04);
	check_count++;
	ma[0] = _mm512_set_epi32(0,2,0,9,0,3,0,9, 0, 2, 0, 9, 0, 3, 0, 9);
	mb[0] = _mm512_set_epi32(0,9,0,9,0,9,0,9, 0, 9, 0, 9, 0, 9, 0, 9);
	b = avx512_neg_cmp_next1(&mask, ma, mb);
	assert(b == true);
	check_count++;

	// check subtract
	__m512i mx[2], my[2], mc;
	mx[0] = _mm512_set_epi32(0,2,0,2,0,3,0,3,0, 2, 0, 2, 0, 3, 0, 3);
	my[0] = _mm512_set_epi32(0,3,0,3,0,2,0,2,0, 3, 0, 3, 0, 2, 0, 2);
	mx[1] = _mm512_set_epi32(0,8,0,9,0,8,0,9,0, 8, 0, 9, 0, 8, 0, 9);
	my[1] = _mm512_set_epi32(0,9,0,8,0,9,0,8,0, 9, 0, 8, 0, 9, 0, 8);
	avx512_subtract2(mx, mx,my);
	e = _mm512_extracti32x8_epi32(mx[0], 0);
	assert(_mm256_extract_epi32(e, 0) == 0x1);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0x1);
	check_count++;
	assert((uint32_t)_mm256_extract_epi32(e, 4) == (uint32_t)0xffffffff);
	check_count++;
	assert((uint32_t)_mm256_extract_epi32(e, 6) == (uint32_t)0xffffffff);
	check_count++;

	// check modular subtract
	mx[0] = _mm512_set_epi32(0,2,0,2,0,3,0,3,0, 2, 0, 2, 0, 3, 0, 3);
	my[0] = _mm512_set_epi32(0,3,0,3,0,2,0,2,0, 3, 0, 3, 0, 2, 0, 2);
	mx[1] = _mm512_set_epi32(0,8,0,9,0,8,0,9,0, 8, 0, 9, 0, 8, 0, 9);
	my[1] = _mm512_set_epi32(0,9,0,8,0,9,0,8,0, 9, 0, 8, 0, 9, 0, 8);
	avx512_modsubtract2(mx,my);
	e = _mm512_extracti32x8_epi32(mx[0], 0);
	assert(_mm256_extract_epi32(e, 0) == 0x1);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0x3);
	check_count++;
	assert((uint32_t)_mm256_extract_epi32(e, 4) == (uint32_t)0xffffffff);
	check_count++;
	assert(_mm256_extract_epi32(e, 6) == 0x2);
	check_count++;

	// check modular subtract with carry
	mx[0] = _mm512_set_epi32(0,2,0,2,0,3,0,3,0, 2, 0, 2, 0, 3, 0, 3);
	my[0] = _mm512_set_epi32(0,3,0,3,0,2,0,2,0, 3, 0, 3, 0, 2, 0, 2);
	mx[1] = _mm512_set_epi32(0,9,0,9,0,9,0,9,0, 9, 0, 9, 0, 9, 0, 9);
	my[1] = _mm512_set_epi32(0,9,0,9,0,9,0,9,0, 9, 0, 9, 0, 9, 0, 9);
	mc    = _mm512_set_epi32(0,0,0,1,0,0,0,1,0, 0, 0, 1, 0, 0, 0, 1);
	avx512_modsubtract32(mx,mc, my);
	e = _mm512_extracti32x8_epi32(mx[0], 0);
	assert(_mm256_extract_epi32(e, 0) == 0x1);
	check_count++;
	assert(_mm256_extract_epi32(e, 2) == 0x1);
	check_count++;
	assert((uint32_t)_mm256_extract_epi32(e, 4) == (uint32_t)0xffffffff);
	check_count++;
	assert(_mm256_extract_epi32(e, 6) == 0x2);
	check_count++;

#endif

	// safeguards 
	b = avxSprpTest(101, 0);
	assert(b == true);
	check_count++;
	b = avxSprpTest(103, 0);
	assert(b == true);
	check_count++;
        b = avxSprpTest(0x000000000000003dull, 0x0000000100000000ull);
	assert(b == true);
	check_count++;


	return check_count;

}

int main(int argc, char **argv)
{
	uint64_t check_count = 0;

	if (argc > 1) {
		printf("%s: sanity checks\n", argv[0]);
		exit(1);
	}

#ifdef __AVX512F__
#else
	printf("AVX512 is NOT enabled\n");
#endif

	printf("Self tests\n");
	check_count += self_tests();

	printf("For each bit size, test pseudoprimes base 2, must be composite\n");
	for (unsigned i = 0; i < sizeof(pseudoprimes) / sizeof(pseudoprimes[0]); i++) {
		uint64_t v_lo;
		uint64_t v_hi;
		bool b;
		uint128_t v = convert128(pseudoprimes[i]);
		v_lo = (uint64_t) v;
		v_hi = (uint64_t) (v >> 64);
		b = avxSprpTest(v_lo, v_hi);
		assert(b == false);
		check_count++;
	}

	printf("For each bit size, test composites base 2, should not escape\n");
	for (unsigned i = 0; i < sizeof(composites) / sizeof(composites[0]); i++) {
		uint64_t v_lo;
		uint64_t v_hi;
		bool b;
		uint128_t v = convert128(composites[i]);
		v_lo = (uint64_t) v;
		v_hi = (uint64_t) (v >> 64);
		b = avxSprpTest(v_lo, v_hi);
		assert(b == false);
		check_count++;
	}

	printf("For each bit size, test non-trivial square roots of 1, should not escape\n");
	for (unsigned i = 0; i < sizeof(square_root) / sizeof(square_root[0]); i++) {
		uint64_t v_lo;
		uint64_t v_hi;
		bool b;
		uint128_t v = convert128(square_root[i]);
		v_lo = (uint64_t) v;
		v_hi = (uint64_t) (v >> 64);
		b = avxSprpTest(v_lo, v_hi);
		assert(b == false);
		check_count++;
	}

	printf("For each bit size, test smallest and largest proven primes, cannot be composite\n");
	for (unsigned i = 0; i < sizeof(primes_min_max) / sizeof(primes_min_max[0]); i++) {
		uint64_t v_lo;
		uint64_t v_hi;
		bool b_min, b_max;
		uint128_t v_min = convert128(primes_min_max[i].p_min);	// smallest prime
		uint128_t v_max = convert128(primes_min_max[i].p_max);	// largest prime
		v_lo = (uint64_t) v_min;
		v_hi = (uint64_t) (v_min >> 64);
		// my_printf(v_min); printf("\n");
		b_min = avx2SprpTest(v_lo, v_hi);
		assert(b_min == true);
		check_count++;
		v_lo = (uint64_t) v_max;
		v_hi = (uint64_t) (v_max >> 64);
		// my_printf(v_max); printf("\n");
		b_max = avxSprpTest(v_lo, v_hi);
		assert(b_max == true);
		check_count++;
	}

	printf("%lu tests passed\n", check_count);
	printf("All tests passed\n");
	return 0;
}