
These tests do functionally check the Fermat and Sprp optimized codes

- sanity checks : from lists of known pseudo-primes, primes, composites, check the correctness of the optimized code output.

- gap_check : display maximal prime gaps found between ranges.

	- this is test code, it only checks for pseudo-maximal pseudo-gaps between pseudo-primes
		This is done to verify that the code under test is MT-safe.
		This is done to verify that the code under test is MT-safe.
	- when OMP is enabled (either at compile time, either at run prime, or both), each thread checks its range separately
		This is done to verify that the code under test is MT-safe.
		I know I should collate the results and send them to a web server on the cloud.
	- progression is about approx 0.3e9 per thread per second a kabylake laptop.
		This is done to verify that the code under test is MT-safe.
		I know I could sieve better, and more, and faster and so on ...
	. gap_check -h for all options

- random_check : checks the algorithms have to correct output, given a bit range to test.

	-start xxx -end yyy [-euler | - sprp] -time zzz

	xxx and yyy are the bit range to test (3 to 128)

	zzz is the test duration in seconds

	-sprp / -euler will be the optimized functions under test.

	- it displays the average processing time every 10 seconds (approx) in rdtsc() cycles or nanoseconds.
