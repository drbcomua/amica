## Install dependencies (Ubuntu)

```shell
sudo apt-get update
sudo apt-get install ocl-icd-opencl-dev g++
```

## Compile

```shell
g++ amicable_gpu.cpp -o amicable_finder -std=c++17 -O3 -pthread -mavx2 -lOpenCL
```

## Run

```shell
# Example: Search for pairs up to 1,000,000
./amicable_finder 1000000

# Example: Search up to 100,000,000 using 64MB segments for the prime sieve
./amicable_finder 100000000 64
```

## Expected result

```shell
$ g++ amicable_gpu.cpp -o amicable_finder -std=c++17 -O3 -pthread -mavx2 -lOpenCL
$ time ./amicable_finder 1000000000 1
Generating primes in memory up to 4294967295 (using 1MB segments, 8 threads)...
Processed segment 256/256
All segments processed. Consolidating primes...
Prime generation complete. Found 203280221 primes.
Using GPU: AMD Radeon R5 M465 Series (radeonsi, iceland, LLVM 19.1.1, DRM 3.57, 6.8.0-60-generic)
Optimizing for GPU: Transferring 3401 primes (up to sqrt(1000000000)) instead of 203280221.
GPU Processing in batches of 2097152 elements.
Processing batch: 998244352 to 1000000000...
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
3,3
948760928=2^5*79*263*1427
951278752=2^5*109*223*1223

4,3
957871508=2^2*11*29*71*97*109
998051692=2^2*11*53*139*3079

Done. Found 562 pairs. Results in amicable_pairs.txt

real    9m3.961s
user    0m46.044s
sys     0m17.773s
```