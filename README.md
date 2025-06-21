# amica

The project dedicated to searching for amicable numbers (https://en.wikipedia.org/wiki/Amicable_numbers).

As task is complex from computational point, and it uses prime numbers, on the first step we'll create utility for
searching prime numbers up to 2^64-1 (which in turn uses results of utility searching prime numbers up to 2^32-1).

## Build amica64 utility

```shell
# build for MacOS, M1 CPU:
g++ -std=c++23 -O3 -ffast-math amica64.cpp -march=armv8.4-a+simd -pthread -o amica64
# build on Windows x86-64 CPU:
g++ -std=c++23 -O3 -ffast-math -pthread -o amica64
```

## Run

```shell
./amica64 100000000
```
This finds all amicable pairs less than 100000000

## Expected result

Tested on Intel Core i5-14400F with 32GB RAM
```shell
$ time ./amica64.exe 1000000000 1
Generating primes in memory up to 4294967295 (using 1MB segments, 16 threads)...
Processed segment 256/256
All segments processed. Consolidating primes...
Prime generation complete. Found 203280221 primes.
Searching for amicable pairs up to 1000000000 using 16 threads
2,1
220=2^2*5*11
284=2^2*71

X2,2
1184=2^5*37
1210=2*5*11^2
...
3,3
948760928=2^5*79*263*1427
951278752=2^5*109*223*1223

4,3
957871508=2^2*11*29*71*97*109
998051692=2^2*11*53*139*3079

Done. Found 564 pairs. Results in amicable_pairs.txt

real    2m44.512s
user    0m0.031s
sys     0m0.031s
```