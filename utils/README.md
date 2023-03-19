# amica/utils
This folder used for prime numbers search utilities and their results.

The following command used to compile uiprimes32:

    g++ -Ofast uiprimes32.cpp -o uiprimes32

The result of running

    ./uiprimes32

is uiprimes32.dat file with all primes less than 2^32.

    Expected size: 813,120,884 bytes
    MD5 sum: 4d052d40fd0ebfa2f015d60062782b0d

Current utility version uses one thread and runs few minutes:
Sample result:

    $ time ./uiprimes32

    real    0m29.350s
    user    0m28.531s
    sys     0m0.793s

64-bit utility
==============

Compile it:

    g++ -Ofast uiprimes64part.cpp -o uiprimes64part

Run as follows:

    ./uiprimes64part 123456

The utility requires uiprimes32.dat file prepared by uiprimes32.
To test that uiprimes64part works properly, run:

    ./uiprimes64part 0
    
The resulting file should be placed at 00/00/00/00.dat and be absolutely the same as uiprimes32.dat (see above).

Expected hash sums:

    4d052d40fd0ebfa2f015d60062782b0d  00/00/00/00.dat
    b35fe0eb9701c802563b5cbc14b42ddd  00/00/00/01.dat
    978a2d15b9c78a45d4cede4a91e41a9d  00/00/00/02.dat
    ...
    56ea92210254169fca00ee8f0bfe1e3a  ff/ff/ff/fe.dat
    a878b16a8990f67abaa80d6e82204c37  ff/ff/ff/ff.dat

Time required to load 32-bit primes, process slice and save findings ranges from 32s to 88s.

Largest prime that can be found using 64-bit utility is 18446744073709551557 (0xffffffffffffffc5)

File sizes:

    00/00/00/00.dat 813120884
    00/00/00/01.dat 761342340
    00/00/00/02.dat 744044304
    ...
    ff/ff/ff/fe.dat 286718248
    ff/ff/ff/ff.dat 228260904
