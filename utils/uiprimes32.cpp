#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#ifdef __GNUC__
#include <immintrin.h>
#endif

// SIMD block: only this function uses AVX2, avoiding ABI mismatches
#ifdef __GNUC__
__attribute__((target("avx2")))
#endif
void collect_primes_simd(const uint64_t* sieve, uint64_t words, uint64_t low, uint64_t seg_bits, std::vector<uint32_t>& outvec) {
    uint64_t w = 0;
#ifdef __AVX2__
    for (; w + 4 <= words; w += 4) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&sieve[w]));
        if (_mm256_testz_si256(v, v)) continue;
        for (int k = 0; k < 4; ++k) {
            uint64_t word = sieve[w + k];
            while (word) {
                int bit = __builtin_ctzll(word);
                uint64_t idx_bit = (w + k) * 64 + bit;
                outvec.push_back(static_cast<uint32_t>(low + 2 * idx_bit));
                word &= word - 1;
            }
        }
    }
#endif
    for (; w < words; ++w) {
        uint64_t word = sieve[w];
        while (word) {
            int bit = __builtin_ctzll(word);
            uint64_t idx_bit = w * 64 + bit;
            outvec.push_back(static_cast<uint32_t>(low + 2 * idx_bit));
            word &= word - 1;
        }
    }
}

int main(int argc, char* argv[]) {
    constexpr uint64_t MAX_N = 0xFFFFFFFFULL;
    // Allow custom segment size in MB via command-line argument (default 32 MB)
    int mb = 32;
    if (argc > 1) {
        mb = std::max(1, std::min(1024, std::atoi(argv[1])));
    }
    uint64_t SEGMENT_BYTES = static_cast<uint64_t>(mb) * 1024 * 1024;
    uint64_t segment_bits = SEGMENT_BYTES * 8;

    // Pre-sieve small primes up to sqrt(MAX_N)
    auto sqrtN = static_cast<uint32_t>(std::sqrt(MAX_N));
    std::vector is_prime_small(sqrtN + 1, true);
    is_prime_small[0] = is_prime_small[1] = false;
    std::vector<uint32_t> primes;
    for (uint32_t i = 2; i <= sqrtN; ++i) {
        if (is_prime_small[i]) {
            primes.push_back(i);
            for (uint64_t j = static_cast<uint64_t>(i) * i; j <= sqrtN; j += i) {
                is_prime_small[j] = false;
            }
        }
    }

    // Calculate total number of odd slots and segments
    uint64_t total_odds = (MAX_N - 3) / 2 + 1;
    auto num_segments = static_cast<uint32_t>((total_odds + segment_bits - 1) / segment_bits);

    // Open output and write the prime '2'
    std::ofstream out("uiprimes32.dat", std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open output file.\n";
        return 1;
    }
    uint32_t two = 2;
    out.write(reinterpret_cast<char*>(&two), sizeof(two));

    // Sequential segmented sieve
    std::vector<uint64_t> sieve;
    std::vector<uint32_t> segment_primes;
    for (uint32_t idx = 0; idx < num_segments; ++idx) {
        uint64_t low_bit = static_cast<uint64_t>(idx) * segment_bits;
        uint64_t high_bit = std::min(low_bit + segment_bits - 1, total_odds - 1);
        uint64_t low = 3 + 2 * low_bit;
        uint64_t seg_bits = high_bit - low_bit + 1;
        uint64_t words = (seg_bits + 63) / 64;

        // Initialize all bits to 'prime'
        sieve.assign(words, ~0ULL);
        if (seg_bits & 63) {
            sieve[words - 1] = ~0ULL >> (64 - (seg_bits & 63));
        }

        // Mark composites
        for (uint32_t p : primes) {
            if (p == 2) continue;
            uint64_t p2 = static_cast<uint64_t>(p) * p;
            if (p2 > low + 2 * high_bit) break;
            uint64_t start = std::max(p2, ( (low + p - 1) / p ) * p);
            if ((start & 1) == 0) start += p;
            uint64_t bit_index = (start - low) / 2;
            for (uint64_t b = bit_index; b < seg_bits; b += p) {
                sieve[b >> 6] &= ~(1ULL << (b & 63));
            }
        }

        // Collect primes in this segment and write them out
        segment_primes.clear();
        segment_primes.reserve(static_cast<size_t>(seg_bits / std::log(low + 1)));
        collect_primes_simd(sieve.data(), words, low, seg_bits, segment_primes);
        out.write(reinterpret_cast<char*>(segment_primes.data()), segment_primes.size() * sizeof(uint32_t));
    }

    out.close();
    std::cout << "Prime generation complete (SIMD, single-threaded, " << mb << "MB segments).\n";
    return 0;
}
