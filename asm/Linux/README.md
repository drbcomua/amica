# amica in assembly language

## Build amica_fast_lnx utility

```shell
# build for Ubuntu, x86-64 CPU:
$ nasm -f elf64 amica_fast_lnx.asm -o amica_fast_lnx.o
$ gcc -no-pie amica_fast_lnx.o -o amica_fast_lnx -lpthread
```

## Run

```shell
./amica_fast_lnx
```
This finds all amicable pairs less than 2^32-1

## Expected result

Tested on Intel Core i5-14400F with 32GB RAM
```shell
$ time ./amica_fast_lnx
Primes found: 203280221, Largest: 4294967291
2,1
220=2^2*5*11
284=2^2*71

X2,2
1184=2^5*37
1210=2*5*11^2

2,2
2620=2^2*5*131
2924=2^2*17*43
...
4,2
4282854730=2*5*7*11*761*7309
5342485430=2*5*1087*491489

4,3
4281566032=2^4*29*97*251*379
4446000368=2^4*113*239*10289

4,4
4280119305=3*5*7*11*59*107*587
4498673655=3*5*7*31*53*89*293

Checked pairs in the range 2..4294967295; found: 1043 amicable pairs

real    2m18.588s
user    36m45.910s
sys     0m1.403s
```
