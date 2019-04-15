# amica/utils
This folder used for prime numbers search utilities and their results.

The following command used to compile uiprimes32:

    g++ -O3 uiprimes32.cpp -o uiprimes32

The result of running

    ./uiprimes32

is uiprimes32.dat file with all primes less than 2^32.

    Expected size: 813,120,884 bytes
    MD5 sum: 4d052d40fd0ebfa2f015d60062782b0d

Current utility version uses one thread and runs few minutes:
Sample result:

    $ time ./uiprimes32

    real    1m9.035s
    user    1m3.031s
    sys     0m4.938s

64-bit utility
==============

Compile it:

    g++ -O3 uiprimes64part.cpp -o uiprimes64part

Run as follows:

    ./uiprimes64part 123456

The utility requires uiprimes32.dat file prepared by uiprimes32.
To test that uiprimes64part works properly, run:

    ./uiprimes64part 0
    
The resulting file should be placed at 00/00/00/00.dat and be absolutely the same as uiprimes32.dat (see above).
