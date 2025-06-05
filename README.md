# amica

The project dedicated to searching for amicable numbers (https://en.wikipedia.org/wiki/Amicable_numbers).

As task is complex from computational point, and it uses prime numbers, on the first step we'll create utility for
searching prime numbers up to 2^64-1 (which in turn uses results of utility searching prime numbers up to 2^32-1).

## Build amica64 utility

```shell
# build for MacOS, M1 CPU:
g++ -std=c++23 -O3 -ffast-math amica64.cpp -march=armv8.4-a+simd -pthread -o amica64
```

## Run

```shell
./amica64 100000000
```
This finds all amicable pairs less than 100000000