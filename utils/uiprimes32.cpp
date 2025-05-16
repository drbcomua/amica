#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>

int main() {
    // Upper bound (2^32 - 1)
    constexpr uint64_t MAX_N = 0xFFFFFFFFULL;
    // Size of each segment in memory (adjust to 32 or 64 MB as needed)
    constexpr uint64_t SEGMENT_BYTES = 32ULL * 1024 * 1024; // 32 MB
    // Each byte here represents an odd number, so range covers twice as many integers
    constexpr uint64_t SEGMENT_RANGE = SEGMENT_BYTES * 2;

    // Pre-sieve up to sqrt(MAX_N)
    uint32_t sqrtN = static_cast<uint32_t>(std::sqrt(MAX_N));
    std::vector<bool> is_prime_small(sqrtN + 1, true);
    is_prime_small[0] = is_prime_small[1] = false;
    std::vector<uint32_t> primes;
    for (uint32_t i = 2; i <= sqrtN; ++i) {
        if (is_prime_small[i]) {
            primes.push_back(i);
            for (uint64_t j = uint64_t(i) * i; j <= sqrtN; j += i) {
                is_prime_small[j] = false;
            }
        }
    }

    // Open output file for writing primes
    std::ofstream out("uiprimes32.dat", std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open output file.\n";
        return 1;
    }

    // Write the first prime (2)
    uint32_t prime = 2;
    out.write(reinterpret_cast<char*>(&prime), sizeof(prime));

    // Segmented sieve for the rest
    for (uint64_t low = 3; low <= MAX_N; low += SEGMENT_RANGE) {
        uint64_t high = std::min(low + SEGMENT_RANGE - 1, MAX_N);
        // Only odd numbers in [low, high]: count = ((high - low) / 2) + 1
        uint64_t segment_size = ((high - low) / 2) + 1;
        std::vector<char> segment(segment_size, 1);

        // Mark non-primes in current segment
        for (uint32_t p : primes) {
            if (p == 2) continue; // skip even prime
            uint64_t p2 = uint64_t(p) * p;
            if (p2 > high) break;
            // Find the first multiple of p within [low, high]
            uint64_t start = (low + p - 1) / p * p;
            if (start < p2) start = p2;
            // Make sure start is odd
            if ((start & 1) == 0) start += p;
            // Mark every 2*p (only odd multiples)
            for (uint64_t j = start; j <= high; j += 2ULL * p) {
                segment[(j - low) / 2] = 0;
            }
        }

        // Write primes from this segment
        for (uint64_t i = 0; i < segment_size; ++i) {
            if (segment[i]) {
                prime = static_cast<uint32_t>(low + 2 * i);
                out.write(reinterpret_cast<char*>(&prime), sizeof(prime));
            }
        }
    }

    out.close();
    std::cout << "Prime generation complete.\n";
    return 0;
}
