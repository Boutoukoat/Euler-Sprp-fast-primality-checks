
# GGG = clang++ -O3 -W -s -m64 -march=native -mno-vzeroupper -flto -static -fomit-frame-pointer -DPARANOID=0 -I . -I ./..
# GGG = aarch64-linux-gnu-g++ -O3 -W -static -DPARANOID=0 -I . -I ./..
# GGG = g++ -O0 -W -g -m64 -march=native -DPARANOID=1 -I . -I ./..
# GGG = clang++ -O0 -W -g -m64 -march=native -DPARANOID=1 -I . -I ./..
# GGG = clang++ -O3 -s -W -m64 -march=native -mavx512f -mavx512dq -mavx512vl -mllvm -x86-use-vzeroupper=0 -flto -static -fomit-frame-pointer -DPARANOID=0 -I . -I ./.. -I ./../tools
GGG = clang++ -O3 -s -W -m64 -march=native -mllvm -x86-use-vzeroupper=0 -flto -static -fomit-frame-pointer -DPARANOID=0 -I . -I ./.. -I ./../tools

all: avx2_gmp_cmp avx2_gmp_perf

avx2_gmp_perf: avx2_gmp_perf.cpp ../m128_utils.h ../m128_utils.cpp ../avx2_sprp.cpp ../avx2_sprp.h ../optimized.h ../tools/divisibility.cpp
	$(GGG) -o avx2_gmp_perf avx2_gmp_perf.cpp ../avx2_sprp.cpp ../m128_utils.cpp ../optimized128.cpp ../optimized65.cpp ../montgomery128.cpp ../tools/divisibility.cpp -lgmp

avx2_gmp_cmp: avx2_gmp_cmp.cpp ../m128_utils.h ../m128_utils.cpp ../avx2_sprp.cpp ../avx2_sprp.h ../optimized.h 
	$(GGG) -o avx2_gmp_cmp avx2_gmp_cmp.cpp ../avx2_sprp.cpp ../m128_utils.cpp ../optimized128.cpp ../optimized65.cpp ../montgomery128.cpp -lgmp

clean:
	rm -f ./avx2_gmp_perf ./avx2_gmp_cmp


