/*
  Utility to calculate all unsigned long primes
*/

#include <iostream>
#include <vector>
#include <inttypes.h>
#include <fstream>

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

inline void save(std::ofstream& f, uint64_t p) {
	uint32_t c = static_cast<uint32_t>(p);
    f.write(reinterpret_cast <char*> (&c), sizeof(uint32_t));
}

int main() {
	std::ofstream fh ("uiprimes32.dat", std::ofstream::out | std::ofstream::binary);
	const uint32_t MAX = 0xFFFFFFFF;
	const uint32_t MAXi = 0xFFFF;
	const uint32_t aSize = 0x100000000/sizeof(uint16_t)/8;
	std::vector<uint16_t> primes(aSize);
	for (uint32_t i=2; i<MAXi; ++i) {
		if(!(primes[i>>4]&mask[i&0b1111])) {
			for (uint64_t j = i*i; j<MAX; j += i) {
				primes[j>>4] |= mask[j&0b1111];
			}
		}
	}
	for (uint64_t i=2; i<MAX; ++i) {
		if(!(primes[i>>4]&mask[i&0b1111])) {
			save(fh, i);
		}
	}
	fh.close();
	return 0;
}