/*
  Utility to calculate primes that exceed uint32_t
  As uint64_t covers extremely large number of primes,
  the calculation takes top 32 bits of the range and writes
  found lower 32 bits of primes into file.
  If AABBCCDD is top 32 bits, then hierarchy is: AA/BB/CC/DD.dat
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <climits>
#include <stdlib.h>
#include <inttypes.h>

using namespace std;

const uint64_t module32 = 0x100000000; // 2^32

inline bool file_is_empty(fstream& f)
{
    return f.peek() == fstream::traits_type::eof();
}

inline void save(fstream& f, uint64_t& c) {
  uint32_t c_mod = c % module32;
  f.write(reinterpret_cast <char*> (&c_mod), sizeof(uint32_t));
}

inline void read(fstream& f, uint64_t& c, const uint32_t& top_bits) {
  uint32_t c_mod;
  f.read(reinterpret_cast <char*> (&c_mod), sizeof(uint32_t));
  c = c_mod + top_bits * module32;
}

int main(int argc, char** argv) {
  if(argc != 2) {
    cout << "Run this utility as follows:\n\t"
         << argv[0]
         << " {slice_number}\nWhere {slice_number} is 0..4294967295 number which is top 32 bits of prime numbers to be found"
         << endl;
    return 1;
  }
  uint32_t slice_number;
  stringstream(argv[1]) >> slice_number;
  short p3 = slice_number % 0x100;
  short p2 = (slice_number >> 8) % 0x100;
  short p1 = (slice_number >> 16) % 0x100;
  short p0 = slice_number >> 24;

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
  uint32_t* prime32 = new uint32_t[num_primes];

  ih.seekg(0);
  ih.read(reinterpret_cast <char*> (prime32), num_primes * sizeof(uint32_t));
  ih.close();

  cout << num_primes << " primes loaded to prime32 array" << endl;

  // then determine largest known uint64_t prime
  if(system(("touch " + fpath.str()).c_str())) {
    cout << "Cannot create slice file:" << fpath.str() << endl;
    return 3;
  }
  fstream oh (fpath.str().c_str(), ios::in | ios::out | ios::binary);
  cout << "Slice file to work with: " << fpath.str().c_str() << endl;

  size_t size64 = 0;

  if(!file_is_empty(oh)) {
    oh.seekg(0, ios::end);
    size64 = oh.tellg();
  }
  uint32_t num_primes_64 = size64/sizeof(uint32_t);

  cout << "64-bit primes file size: " << size64 << "; number of stored primes: " << num_primes_64 << endl;

  uint64_t last_prime;
  if(num_primes_64==0) {
    if(slice_number) {
      last_prime = slice_number * module32 - slice_number * module32 % 6 + 1; // first odd before slice
    } else { // 0 slice should start from 5 as 2 and 3 are not in 6k+-1 sieve
      oh.clear();
      oh.seekp(0, ios::beg); // put cursor to the beginning of empty file
      save(oh, last_prime = 2);
      save(oh, last_prime = 3);
      save(oh, last_prime = 5);
    }
    cout << "No 64-bit primes known yet, so setting last_prime=" << last_prime << endl;
  } else {
    oh.seekg((num_primes_64 - 1) * sizeof(uint32_t), ios::beg);
    read(oh, last_prime, slice_number);
    cout << "Largest 64-bit prime known: last_prime=" << last_prime << endl;
  }
  oh.clear();
  oh.seekp(0, ios::end); // put cursor to the end

  uint64_t k = (last_prime + 1) / 6;
  short l = last_prime - 6 * k;
  
  uint64_t candidate;

  cout << "Last known prime: " << last_prime << " = 6*" << k << (l<0?"-":"+") << "1" << endl;

  lbl_1:
  while(1) {
    if ((l = -l) == -1) ++k;
    candidate = 6 * k + l; // candidate prime

    if (slice_number != candidate / module32) break; // job's done for this slice

    size_t i = 0; // prime index
    uint64_t p; // test known prime as divisor
    do {
      if(!(candidate % (p = prime32[++i]))) goto lbl_1;
    } while ((p * p < candidate) && (i < num_primes - 1));

    i = 0; // at this point i counts 64-bit primes

    save(oh, candidate);
    oh.flush();
  }

  oh.close();
  cout << "Search finished" << endl;

  return 0;
}