

Optimized code for Fermat base 2 and SPRP base 2 primality checks for odd numbers n from 3 to 128 bits.
	n >= 0x3
	n < 0x100000000000000000000000000000000

Fermat tests includes the Euler criterion and a fast way to compute it (code is branch-free).
https://en.wikipedia.org/wiki/Fermat_primality_test
https://en.wikipedia.org/wiki/Euler_pseudoprime
https://en.wikipedia.org/wiki/Euler%27s_criterion
https://en.wikipedia.org/wiki/Legendre_symbol

Sprp tests are based on this code
https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test

Preferred entry points for all values of n_hi and n_lo

	bool optimizedSprpTest(uint64_t n_lo, uint64_t n_hi);
	bool optimizedFermatTest(uint64_t n_lo, uint64_t n_hi);

The required files are:
	montgomery.h optimized.h m128_utils.h m64_utils.h m_reg.h
	montgomery128.cpp  optimized128.cpp  optimized65.cpp

Example of code to show how to use the tests:
	tests/gap_check.cpp

SPRP tests are, by design, not really branch-free, but have less pseudoprimes than fermat test.
This has to be taken into account on targets where conditional unpredictable branchs are expensive at run time.

The following files implement different variants of the same algorithm.

- montgomery128.cpp : montgomery CIOS textbook version

	https://thomas-plantard.github.io/pdf/Plantard21.pdf
	https://eprint.iacr.org/2016/487.pdf   (algorithm 2)

- optimized128.cpp : optimized code, 

	less textbook than it should be.

- optimized65.cpp : 65-bits primality check

	Optimized down to the very last cycle for odd numbers like 0x1xxxxxxxxxxxxxxxx
	Source code tend to be more complex to read and understand.

- tools/slow.cpp : very slow C code, provably correct against overflows.

	this code is for debug purposes. It is the reference test code.

- tools/generic.cpp : reasonably slow C code 

	suitable for test purposes

- tools/barrett.cpp : reasonably fast C code 

	suitable for test purposes

- Makefile

	This is an example to show how to build the code on an unix-like systems.

		$ make all

	There is no library, since library building is non portable across compilers and
	architectures. Please feel free to build your own library the way you prefer it 
	(dynamic, static, 1 single capsule file, object files, w/ assert included or excluded 
	.lib / .a / .o / .obj / .so / .dll / .exe / .bin ...)

- Directories for test code

	./tests
	./unit_tests
	./perf

	To build and quicky test the code and the test code (stops on first failure)
		$ make clean
		$ make check

	To run the perf tests and get csv files for Excel, stored in ./misc directory
		$ make csv

This code for gcc and clang compilers makes heavy use of x86-64 assembly language. 
To use it on a different compiler and/or a different hw architecture, toggle INLINE_ASM from 1 to 0 in m_reg.h

__arch64__ support:

	I don't have an ARM under my heavy hands, I did install an emulation for Debian x86-64

	sudo apt-get install qemu binfmt-support qemu-user-static qemu-user qemu-user-static \
 		gcc-aarch64-linux-gnu g++-aarch64-linux-gnu \
		binutils-aarch64-linux-gnu binutils-aarch64-linux-gnu-dbg gdc-aarch64-linux-gnu \
		build-essential

	Then, all is needed is to enable this line in Makefile

	GGG = aarch64-linux-gnu-g++ -O3 -W -static -DPARANOID=1 -I . -I ./tools

__avx2__ support:
__avx512__ support:

	avx2 is used to run 4 sprp tests at the same time for deterministic primality checks
	avx512 is used to run 8 sprp tests at the same time for deterministic primality checks

__ifma__ support:
	(work in progress) avx512-ifma is used to run 4 sprp tests at the same time for deterministic primality checks

uint128_t support:
Note that the compiler must support natively a uint128_t type, and such a class implementation 
can be found on the web, like https://github.com/calccrypto/uint128_t (first google result). 
Note also that LLVM supports uint128_t since 4.1 and clang 3.1 since about 2006. Nothing really new ....

C23, .NET 7 have their own definition of the thing. 
openCl seems to implement 128 bits with unsigned long long.
VS2013 and VS2015 documentation mentions __int128 and some MSVC header files might have std::_Signed128 or std::_Unsigned128.


-----------------------------------------------------------------------------------------------
LICENSE
-----------------------------------------------------------------------------------------------

The MIT License (MIT)

Copyright (c) 2024

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

