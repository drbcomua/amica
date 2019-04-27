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
#include <vector>
#include <stdlib.h>
#include <inttypes.h>

using namespace std;

// Mask to translate lower 4 bits to 'hot unit'
const uint16_t mask[] = {
	0b0000000000000001,
	0b0000000000000010,
	0b0000000000000100,
	0b0000000000001000,
	0b0000000000010000,
	0b0000000000100000,
	0b0000000001000000,
	0b0000000010000000,
	0b0000000100000000,
	0b0000001000000000,
	0b0000010000000000,
	0b0000100000000000,
	0b0001000000000000,
	0b0010000000000000,
	0b0100000000000000,
	0b1000000000000000
};

// Mask to initialize array of primes (as we know 2 is prime)
const uint16_t mask2 = 0b0101010101010101;

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
  fstream oh (fpath.str().c_str(), ios::out | ios::binary);
  cout << "Slice file to work with: " << fpath.str().c_str() << endl;

    const uint32_t aSize = 0x100000000/sizeof(uint16_t)/8; // size of 32-bit range packed into 16-bit words
	std::vector<uint16_t> primes(aSize);

	// Init by excluding even numbers
	for (uint32_t i=0; i<aSize; ++i) {
		primes[i] = mask2;
	}
  if (!slice_number) {
    primes[0] &= 0b1111111111111011;
    primes[0] |= 0b11;
  }

	uint32_t i=1;
	do {
      uint64_t p = prime32[i];
			uint64_t j = slice_number ? p - slice_number*module32%p : p*p;
      if(!(j%2)) j+=p;
			while (j <= UINT32_MAX) {
				primes[j>>4] |= mask[j&0b1111];
        j += p+p;
			}
	} while (++i < num_primes);

	for (uint64_t i=0; i<=UINT32_MAX; ++i) { // save found primes
		if(!(primes[i>>4]&mask[i&0b1111])) {
			uint32_t c = static_cast<uint32_t>(i);
			oh.write(reinterpret_cast <char*> (&c), sizeof(uint32_t));
		}
	}
	oh.close();

  return 0;
}
