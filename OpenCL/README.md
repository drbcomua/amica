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