# amica/utils
This folder used for prime numbers search utilities and their results.

The following command used to compile uiprimes32 on x86-64:

    g++ -std=c++23 -Ofast -mavx2 uiprimes32.cpp -o uiprimes32

On Apple M1:

     g++ -std=c++23 -O3 -ffast-math uiprimes32.cpp -march=armv8.4-a+simd -o uiprimes32

The result of running

    ./uiprimes32

is uiprimes32.dat file with all primes less than 2^32.

    Expected size: 813,120,884 bytes
    MD5 sum: 4d052d40fd0ebfa2f015d60062782b0d

Current utility version uses one thread and runs few minutes:
Sample result:

    % time ./uiprimes32 4
    Prime generation complete (SIMD, single-threaded, 4MB segments).
    ./uiprimes32 4  7.04s user 0.44s system 98% cpu 7.561 total

where 4 is segment size (4MB) - optimal size differ for different CPUs.

64-bit utility
==============

Compile it on x86-64:

    g++ -std=c++23 -Ofast -mavx2 uiprimes64part.cpp -o uiprimes64part

On Apple M1:

    g++ -std=c++23 -O3 -ffast-math uiprimes64part.cpp -march=armv8.4-a+simd -o uiprimes64

Run as follows:

    ./uiprimes64 123456

The utility requires uiprimes32.dat file prepared by uiprimes32.
To test that uiprimes64 works properly, run:

    ./uiprimes64 0
    
The resulting file should be placed at 00/00/00/00.dat and be absolutely the same as uiprimes32.dat (see above).

Expected hash sums:

    4d052d40fd0ebfa2f015d60062782b0d  00/00/00/00.dat
    b35fe0eb9701c802563b5cbc14b42ddd  00/00/00/01.dat
    978a2d15b9c78a45d4cede4a91e41a9d  00/00/00/02.dat
    ...
    55cf2684a7b90fdf6172ac8fa2135d86  FF/FF/FF/FE.dat
    c25050b5aefb705a2b89df0c02400b39  FF/FF/FF/FF.dat

Time required to load 32-bit primes, process slice and save findings ranges from 32s to 88s.

Largest prime that can be found using 64-bit utility is 18446744073709551557 (0xffffffffffffffc5)

File sizes:

    00/00/00/00.dat 813120884
    00/00/00/01.dat 761342340
    00/00/00/02.dat 744044304
    ...
    FF/FF/FF/FE.dat 387254980
    FF/FF/FF/FF.dat 387192372


Verify utility
==============

Utility runs deterministic Miller–Rabin primality test.
By default it checks all found numbers in FF/FF/FF/FF.dat file.