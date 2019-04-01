# amica/utils
This folder used for prime numbers search utilities and their results.

The following command used to compile uiprimes32:

    g++ -O3 uiprimes32.cpp -o uiprimes32

The result of running

    ./uiprimes32

is uiprimes32.dat file with all primes less than 2^32.

    Expected size: 813,120,884 bytes
    MD5 sum: 4d052d40fd0ebfa2f015d60062782b0d

Current utility version uses one thread and runs 2..4 hours depending on processor core speed.
Sample result:

    $ time ./uiprimes32

    real    141m46.326s
    user    140m42.063s
    sys     0m19.453s

64-bit utility
==============

Compile it:

    g++ -O3 uiprimes64part.cpp -o uiprimes64part

Run as follows:

    ./uiprimes64part 123456
