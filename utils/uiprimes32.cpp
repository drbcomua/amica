/*
  Utility to calculate all unsigned long primes
*/

#include <iostream>
#include <vector>
#include <inttypes.h>
#include <fstream>

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

int main() {
	std::ofstream fh ("uiprimes32.dat", std::ofstream::out | std::ofstream::binary);
	const uint32_t aSize = 0x100000000/sizeof(uint16_t)/8; // size of 32-bit range packed into 16-bit words
	std::vector<uint16_t> primes(aSize);
	
	// Init by exluding even numbers except 2
	for (uint32_t i=0; i<aSize; ++i) {
		primes[i] = mask2;
	}
	primes[0] &= 0b1111111111111011;
	
	uint32_t i=3;
	do { // Check divisor up to sqrt(UINT32_MAX) = UINT16_MAX
		if(!(primes[i>>4]&mask[i&0b1111])) {
			// bybass even numbers and start from squared prime
			uint64_t j = i*i;
			do {
				primes[j>>4] |= mask[j&0b1111];
			} while ((j += i+i) < UINT32_MAX);
		}
	} while ((i += 2) < UINT16_MAX);
	
	for (uint64_t i=2; i<UINT32_MAX; ++i) { // save found primes
		if(!(primes[i>>4]&mask[i&0b1111])) {
			uint32_t c = static_cast<uint32_t>(i);
			fh.write(reinterpret_cast <char*> (&c), sizeof(uint32_t));
		}
	}
	fh.close();
	return 0;
}