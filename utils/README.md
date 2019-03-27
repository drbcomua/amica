# amica/utils
This folder used for prime numbers search utilities and their results.

The following command used to compile uiprimes32:
    g++ -O6 uiprimes32.cpp -o uiprimes32

The result of running
    ./uiprimes32
is uiprimes32.dat file with all primes less than 2^32.

Current utility version uses one thread and runs 2..4 hours depending on processor core speed.
Sample result:
    $ time ./uiprimes32

    real    141m46.326s
    user    140m42.063s
    sys     0m19.453s