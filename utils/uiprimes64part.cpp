/*
  Utility to calculate primes that exceed uint32_t
  As uint64_t covers extremely large number of primes,
  the calculation takes top 32 bits of the range and writes
  found lower 32 bits of primes into file.
  If AABBCCDD is top 32 bits, then hierarchy is: AA/BB/CC/DD.dat
*/

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

inline uint32_t index2int(uint32_t i) {
    return (i<<1)+1;
}

inline uint32_t int2index(uint32_t i) {
    return (i-1)>>1;
}

inline void write(ofstream &fh, uint32_t p) {
    fh.write(reinterpret_cast <char*> (&p), sizeof(uint32_t));
}

int main(int argc, char** argv) {
    const uint32_t MAX_VALUE = 0xFFFFFFFF;
    const uint32_t MAX_INDEX = int2index(MAX_VALUE);

    vector<bool> primes(MAX_INDEX, true);

    if(argc != 2) {
        cout << "Run this utility as follows:\n\t"
             << argv[0]
             << " {slice_number}\nWhere {slice_number} is 0..4294967295 number which is top 32 bits of prime numbers to be found"
             << endl;
        return 1;
    }
    uint64_t slice_number;
    stringstream(argv[1]) >> slice_number;
    unsigned short p3 = slice_number % 0x100;
    unsigned short p2 = (slice_number >> 8) % 0x100;
    unsigned short p1 = (slice_number >> 16) % 0x100;
    unsigned short p0 = slice_number >> 24;

    ostringstream path;
    path << hex
         << setfill('0') << setw(2) << p0 << "/"
         << setfill('0') << setw(2) << p1 << "/"
         << setfill('0') << setw(2) << p2 << "/";
    cout << path.str() << endl;

    ostringstream fpath;
    fpath << path.str() << hex << setfill('0') << setw(2) << p3 << ".dat";

    if(system(("mkdir -p " + path.str()).c_str())) {
        cout << "Cannot create folder: " << path.str() << endl;
        return 2;
    }

    // first load known uint32_t primes
    cout << "Loading primes..." << endl;
    ifstream ih ("uiprimes32.dat", ios::in | ios::binary);

    ih.seekg(0, ios::end);
    size_t size = ih.tellg();
    uint32_t num_primes = size/sizeof(uint32_t) - 1;
    auto* prime32 = new uint32_t[num_primes];

    ih.seekg(0);
    ih.read(reinterpret_cast <char*> (prime32), num_primes * sizeof(uint32_t));
    ih.close();

    cout << num_primes << " primes loaded to prime32 array" << endl;

    if(system(("touch " + fpath.str()).c_str())) {
        cout << "Cannot create slice file:" << fpath.str() << endl;
        return 3;
    }

    ofstream oh (fpath.str().c_str(), ofstream::out | ofstream::binary);
    cout << "Slice file to work with: " << fpath.str().c_str() << endl;

    for (uint32_t i=1; i < num_primes; ++i) { // skipping 2
        uint64_t d = prime32[i]; // prime divisor
        if ((d*d)>>32 > slice_number) {
            break; // divisor is too big for current slice
        }
        uint64_t mod = (slice_number<<32)%d;
        uint32_t shift = slice_number?int2index(d-mod+d*(mod%2)):int2index(d)+d;
        for (uint32_t j=shift; j<=MAX_INDEX; j+=d) {
            primes[j]=false;
        }
    }

    if (slice_number == 0) {
        write(oh, 2);
        primes[0]=false; // 1 is not prime
    }
    for (uint64_t i=0; i<MAX_INDEX; ++i) {
        if (primes[i]) {
            write(oh, index2int(i));
        }
    }

    oh.close();
    return 0;
}