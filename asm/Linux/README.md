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

real    0m8.473s
user    2m14.254s
sys     0m0.230s
```

## Build amica_fast_lnx_40 utility

Same algorithm extended to the 2^40 range. Needs ~540 MB of RAM.

```shell
# build for Ubuntu, x86-64 CPU:
$ nasm -f elf64 amica_fast_lnx_40.asm -o amica_fast_lnx_40.o
$ gcc -no-pie amica_fast_lnx_40.o -o amica_fast_lnx_40 -lpthread
```

## Run

```shell
./amica_fast_lnx_40 -n
```
This finds all amicable pairs less than 2^40. The optional `-n` flag skips
the prime-counting wave, which costs about a minute at this range and
contributes nothing to the pair search; drop it to get the
`Primes found:` line as well.

## Expected result

Tested on Intel Core i5-14400F with 32GB RAM
```shell
$ time ./amica_fast_lnx_40 -n
...
5,4
1095964542310=2*5*11*53*79*167*14249
1137924737690=2*5*13*89*683*143999

4,3
1097623165180=2^2*5*17*61*227*233141
1393922490692=2^2*61*293*683*28547

3,3
1097676143685=3^2*5*11*19*151*359*2153
1108791849915=3^2*5*17*19*37*359*5743

Checked pairs in the range 2..1099511627775; found: 7893 amicable pairs

real    89m35.448s
user    1430m34.357s
sys     0m16.508s
```
