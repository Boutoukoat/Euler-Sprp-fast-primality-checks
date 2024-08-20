
These tests do run some simplistic performance tests of some of the algorithms, together with a check of the result correctness.

fermat_perf.cpp  : 
	- fastest algorithms with fermat base 2 test + euler criterion
sprp_perf.cpp :
	- fastest algorithms with sprp base 2 test

generic_perf.cpp :
	- textbook algorithms with fermat test
		- barrett : code based on barrett reduction
		- generic : code based on << and %. This is quite slow, and most of the code is easy to prove.
		- slow : simple generic code based on modular add. This is badly slow and all the code is provable.


for these programs, the command line options are

	-start xxx -end yyy [-slow] [-csv | -nocsv]

where xxx is as low as 3 and yyy is as high as 128, to check the bit-range intervals and display the performance, in rdtsc() wall-clock 
cycles, of the different algorithms.

-slow option adds the slowest and provably unbreakable variant where the only operation is a modular addition. 
No overflow or glitch is possibly possible with that code. Proven correctness comes at a serious cycle cost.

-csv option enables a csv output with can be saved in a file compatible with Excel input.




