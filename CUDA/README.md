## Install dependencies

Download and install CUDA Toolkit

## Compile

```shell
nvcc amicable_cuda.cu -o amicable_cuda -std=c++17 -O3 -Xcompiler "-pthread"
```

## Run

```shell
# Example: Search for pairs up to 1,000,000
./amicable_cuda 1000000

# Example: Search up to 100,000,000 using 64MB segments for the prime sieve
./amicable_cuda 100000000 64
```

## Expected result

```shell
$ time ./amicable_cuda.exe 1000000000 1
Generating primes in memory up to 4294967295 (using 1MB segments, 16 threads)...
Processed segment 256/256
All segments processed. Consolidating primes...
Prime generation complete. Found 203280221 primes.
Using GPU: NVIDIA GeForce RTX 4070 SUPER
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

Done. Found 564 pairs. Results in amicable_pairs.txt

real    0m19.548s
user    0m0.030s
sys     0m0.000s
```