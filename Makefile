# GGG = aarch64-linux-gnu-g++ -O3 -W -static -DPARANOID=0 -I . -I ./tools
# GGG = clang++ --target=aarch64-pc-linux -O3 -W -static -DPARANOID=0 -I . -I ./tools   # untested !
# GGG = g++ -O0 -W -g -m64 -march=native -DPARANOID=1 -I . -I ./tools
# GGG = clang++ -O0 -W -g -m64 -march=native -DPARANOID=1 -I . -I ./tools
# GGG = g++ -O3 -W -s -m64 -march=native -fomit-frame-pointer -fexpensive-optimizations -DPARANOID=0 -I . -I ./tools
GGG = clang++ -O3 -W -s -m64 -march=native -flto -fomit-frame-pointer -DPARANOID=0 -I . -I ./tools

all: all_unit_tests all_tests all_perf

all_unit_tests: unit_tests/m_reg_tests unit_tests/generic_tests unit_tests/barrett_tests unit_tests/optimized65_tests

all_tests: tests/sanity_check tests/random_test tests/gap_check

all_perf: perf/fermat_perf perf/sprp_perf perf/generic_perf

perf/sprp_perf: perf/sprp_perf.cpp m128_utils.h m128_utils.cpp m_reg.h optimized128.cpp optimized65.cpp montgomery128.cpp montgomery.h tools/slow.cpp tools/slow.h tools/generic.cpp tools/generic.h
	$(GGG) -o perf/sprp_perf perf/sprp_perf.cpp m128_utils.cpp optimized128.cpp optimized65.cpp montgomery128.cpp tools/slow.cpp tools/generic.cpp

perf/generic_perf: perf/generic_perf.cpp m128_utils.h m128_utils.cpp m_reg.h tools/slow.cpp tools/slow.h tools/generic.cpp tools/generic.h tools/barrett.cpp tools/barrett.h
	$(GGG) -o perf/generic_perf perf/generic_perf.cpp m128_utils.cpp tools/slow.cpp tools/generic.cpp tools/barrett.cpp

perf/fermat_perf: perf/fermat_perf.cpp m128_utils.h m128_utils.cpp m_reg.h optimized128.cpp optimized65.cpp montgomery128.cpp montgomery.h tools/slow.cpp tools/slow.h tools/generic.cpp tools/generic.h
	$(GGG) -o perf/fermat_perf perf/fermat_perf.cpp m128_utils.cpp optimized128.cpp optimized65.cpp montgomery128.cpp tools/slow.cpp tools/generic.cpp

tests/sanity_check: tests/sanity_check.cpp m128_utils.h m128_utils.cpp m_reg.h optimized128.cpp optimized65.cpp montgomery128.cpp montgomery.h
	$(GGG) -o tests/sanity_check tests/sanity_check.cpp m128_utils.cpp optimized128.cpp optimized65.cpp montgomery128.cpp

tests/gap_check: tests/gap_check.cpp m128_utils.h m128_utils.cpp m_reg.h optimized128.cpp optimized65.cpp montgomery128.cpp montgomery.h tests/gap_check.wheel
	$(GGG) -fopenmp=libomp -o tests/gap_check tests/gap_check.cpp m128_utils.cpp optimized128.cpp optimized65.cpp montgomery128.cpp

tests/random_test: tests/random_test.cpp m128_utils.h m128_utils.cpp m_reg.h optimized128.cpp optimized65.cpp montgomery128.cpp montgomery.h tools/generic.cpp tools/generic.h tools/slow.cpp tools/slow.h
	$(GGG) -o tests/random_test tests/random_test.cpp m128_utils.cpp optimized65.cpp optimized128.cpp montgomery128.cpp tools/generic.cpp tools/slow.cpp

unit_tests/m_reg_tests: unit_tests/m_reg_tests.cpp m_reg.h
	$(GGG) -o unit_tests/m_reg_tests unit_tests/m_reg_tests.cpp

unit_tests/generic_tests: unit_tests/generic_tests.cpp tools/generic.cpp m128_utils.h m128_utils.cpp tools/barrett.cpp tools/barrett.h tools/slow.cpp tools/slow.h
	$(GGG) -o unit_tests/generic_tests unit_tests/generic_tests.cpp m128_utils.cpp tools/barrett.cpp tools/slow.cpp

unit_tests/barrett_tests: unit_tests/barrett_tests.cpp tools/generic.cpp m128_utils.h m128_utils.cpp tools/barrett.cpp tools/barrett.h tools/slow.cpp tools/slow.h
	$(GGG) -o unit_tests/barrett_tests unit_tests/barrett_tests.cpp m128_utils.cpp tools/slow.cpp

unit_tests/optimized65_tests: unit_tests/optimized65_tests.cpp optimized65.cpp m128_utils.h m128_utils.cpp montgomery128.cpp
	$(GGG) -o unit_tests/optimized65_tests unit_tests/optimized65_tests.cpp m128_utils.cpp

clean:
	rm -f ./perf/fermat_perf ./perf/sprp_perf ./perf/generic_perf 
	rm -f ./tests/random_test ./tests/sanity_check ./tests/gap_check
	rm -f ./unit_tests/generic_tests unit_tests/m_reg_tests ./unit_tests/barrett_tests ./unit_tests/optimized65_tests

csv: perf/fermat_perf perf/sprp_perf perf/generic_perf 
	perf/fermat_perf -csv | tee misc/fermat.csv
	perf/sprp_perf -csv | tee misc/sprp.csv
	perf/generic_perf -csv | tee misc/generic.csv

zip:
	zip -r pseudoprimes.zip \
		unit_tests/ tests/ perf/ misc/ tools/ \
		*.cpp *.h Readme.txt Makefile \
		unit_tests/*.cpp unit_tests/Readme.txt \
		tests/*.cpp tests/Readme.txt tests/*.wheel \
		tools/*.cpp tools/*.h tools/Readme.txt \
		perf/*.cpp perf/Readme.txt

check: tests/sanity_check tests/random_test \
	perf/fermat_perf perf/sprp_perf perf/generic_perf \
	unit_tests/m_reg_tests unit_tests/generic_tests unit_tests/barrett_tests unit_tests/optimized65_tests
	# ----------------------------------------------------------
	# Unit tests
	# ----------------------------------------------------------
	./unit_tests/m_reg_tests
	./unit_tests/generic_tests
	./unit_tests/barrett_tests
	./unit_tests/optimized65_tests
	# ----------------------------------------------------------
	# Sanity checks
	# ----------------------------------------------------------
	./tests/sanity_check
	# ----------------------------------------------------------
	# Fermat-euler tests
	# ----------------------------------------------------------
	./perf/fermat_perf -start 63 -end 70 
	# ----------------------------------------------------------
	# Sprp tests
	# ----------------------------------------------------------
	./perf/sprp_perf -start 63 -end 70
	# ----------------------------------------------------------
	# Random Fermat tests for 60 seconds
	# ----------------------------------------------------------
	./tests/random_test -start 3 -end 128 -time 60
	# ----------------------------------------------------------
	# Random Fermat tests for 20 seconds
	# ----------------------------------------------------------
	./tests/random_test -start 65 -end 65 -time 20
	# ----------------------------------------------------------
	@echo "That's all, folks !"

