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
% time ./amicable_metal 100000000 2 
Generating primes in memory up to 4294967295 (using 2MB segments, 8 threads)...
Processed segment 128/128
All segments processed. Consolidating primes...
Prime generation complete. Found 203280221 primes.
Using GPU: Apple M1
Optimizing for Metal: Using 1229 primes (up to sqrt(100000000)) in shared memory.
GPU Processing in batches of 2097152 elements.
Processing batch: 98566144 to 100000000...
All batches processed.
Filtering potential pairs on CPU...
GPU computation and filtering complete.
2,1
220=2^2*5*11
284=2^2*71

X2,2
1184=2^5*37
1210=2*5*11^2
...
2,2
95629904=2^4*37*67*2411
97580944=2^4*67*227*401

2,3
96304845=3*5*7^2*13*10079
96747315=3*5*7^2*23*59*97

2,2
97041735=3^2*5*7*71*4339
97945785=3^2*5*7*239*1301

Done. Found 231 pairs. Results in amicable_pairs.txt
./amicable_metal 100000000 2  9.55s user 0.83s system 16% cpu 1:02.13 total
```