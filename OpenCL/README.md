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
# Example: Search for pairs with lower number up to 1,000,000
./amicable_finder 1000000
```

## Expected result

```shell
$ g++ amicable_gpu.cpp -o amicable_finder -std=c++17 -O3 -pthread -mavx2 -lOpenCL
$ time ./amicable_finder 1000000000
Generating primes up to 31623...
Prime generation complete. Found 3401 primes.
Using OpenCL device: AMD Radeon R5 M465 Series (radeonsi, iceland, LLVM 19.1.1, DRM 3.57, 6.8.0-62-generic)
Building OpenCL kernel...
Transferring 3401 primes to the device.
Processing numbers up to 1000000000 on the device...
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
X2,2
992776995=3^2*5*7^2*450239
1008990045=3^2*5*7*797*4019

X4,2
994945490=2*5*7^2*11^2*97*173
1331936326=2*7^2*881*15427

4,5
996088412=2^2*11*17*151*8819
1030959268=2^2*23*37*41*83*89

Done. Found a total of 586 pairs. Results also saved in amicable_pairs.txt

real	18m0.280s
user	9m15.848s
sys	0m19.229s
```