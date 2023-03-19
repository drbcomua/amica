/*
  Utility to calculate all unsigned long primes
*/

#include <vector>
#include <fstream>

using namespace std;

inline uint32_t index2int(uint32_t i) {
    return (i << 1) + 3;
}

inline uint32_t int2index(uint32_t i) {
    return (i - 3) >> 1;
}

inline void write(ofstream &fh, uint32_t p) {
    fh.write(reinterpret_cast <char *> (&p), sizeof(uint32_t));
}

int main() {
    const uint32_t MAX_VALUE = 0xFFFFFFFF;
    const uint32_t MAX_INDEX = int2index(MAX_VALUE);
    const uint32_t SQRT_MAX_INDEX = int2index(MAX_VALUE >> 16);

    std::vector<bool> primes(MAX_INDEX, true);

    ofstream fh("uiprimes32.dat", ofstream::out | ofstream::binary);
    write(fh, 2);

    for (uint32_t i = 0; i < SQRT_MAX_INDEX; ++i) {
        if (primes[i]) {
            uint32_t prime = index2int(i);
            write(fh, prime);
            for (uint64_t j = prime * prime; j <= MAX_VALUE; j += (prime << 1)) {
                primes[int2index(j)] = false;
            }
        }
    }

    for (uint32_t i = SQRT_MAX_INDEX; i < MAX_INDEX; ++i) {
        if (primes[i]) {
            write(fh, index2int(i));
        }
    }

    fh.close();
    return 0;
}