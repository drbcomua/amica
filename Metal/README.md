## Install dependencies (MacOS)

Download metal-cpp from https://developer.apple.com/metal/cpp/ (should correspond to OS version)

## Compile

```shell
clang++ amicable_metal.cpp -o amicable_metal \
        -std=c++17 -O3 -pthread \
        -I./metal-cpp \
        -framework Foundation -framework Metal -framework QuartzCore
```

## Run

```shell
# Example: Search for pairs up to 10,000,000
./amicable_metal 10000000
```

## Expected result

```shell
% clang++ -std=c++17 -O3 -o amicable_metal amicable_metal.cpp -I./metal-cpp -framework Foundation -framework Metal -framework QuartzCore
% time ./amicable_metal 1000000000
Generating primes up to 31623...
Prime generation complete. Found 3401 primes.
Using GPU: Apple M1
Transferring 3401 primes to the device.
GPU Processing in batches of 8388608 elements.
Processing batch: 998244353 to 1000000000...
All batches processed.

--- Found 586 Amicable Pairs ---

2,1
220=2^2*5*11
284=2^2*71

X2,2
1184=2^5*37
1210=2*5*11^2
...
X4,2
994945490=2*5*7^2*11^2*97*173
1331936326=2*7^2*881*15427

4,5
996088412=2^2*11*17*151*8819
1030959268=2^2*23*37*41*83*89

Done. Found a total of 586 pairs. Results also saved in amicable_pairs.txt
./amicable_metal 1000000000  74.57s user 2.59s system 4% cpu 28:44.76 total
```