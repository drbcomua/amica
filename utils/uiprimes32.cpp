/*
  Utility to calculate all unsigned long primes
*/

#include <iostream>
#include <fstream>
#include <climits>
#include <vector>

using namespace std;

inline void save(ofstream& f, unsigned& c) {
  f.write(reinterpret_cast <char*> (&c), sizeof(unsigned));
}

int main() {
  ofstream fh ("uiprimes32.dat", ofstream::out | ofstream::binary);
// To Do: check if uint32_t works better on different platforms
  unsigned k = 0;
  short l = 1;
  const unsigned MAXK = (UINT_MAX - 1) / 6;
  unsigned candidate;

  vector<unsigned> primes;
  primes.push_back(candidate = 2);
  save(fh, candidate);
  primes.push_back(candidate = 3);
  save(fh, candidate);

  lbl_1:
  while(k < MAXK) {
    if ((l = -l) == -1) ++k;
    unsigned candidate = 6 * k + l; // candidate prime

    size_t i = 0;    // prime index
    unsigned long p; // test known prime as divisor
    do {
      if(!(candidate % (p = primes[++i]))) goto lbl_1;
    } while (p * p < candidate);

    primes.push_back(candidate); //no known divisors found, so this is new prime
    save(fh, candidate);
  }

  fh.close();

  return 0;
}