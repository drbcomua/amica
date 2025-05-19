#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#if defined(__x86_64__)
#  include <immintrin.h>    // AVX2
#endif
#if defined(__aarch64__)
#  include <arm_neon.h>     // NEON
#endif

#if defined(__x86_64__)
__attribute__((target("avx2")))
#endif
void collect_primes_simd(const uint64_t* sieve,
                         uint64_t words,
                         uint64_t low,
                         uint64_t seg_bits,
                         std::vector<uint32_t>& out_vec)
{
    uint64_t w = 0;
#ifdef __AVX2__
    for (; w + 4 <= words; w += 4) {
        __m256i v = _mm256_loadu_si256((__m256i*)&sieve[w]);
        if (_mm256_testz_si256(v, v)) continue;
        for (int k = 0; k < 4; ++k) {
            uint64_t word = sieve[w + k];
            while (word) {
                int bit = __builtin_ctzll(word);
                uint64_t idx = (w + k)*64 + bit;
                out_vec.push_back(uint32_t(low + 2*idx));
                word &= word - 1;
            }
        }
    }
#elif defined(__aarch64__)
    for (; w + 2 <= words; w += 2) {
        uint64x2_t v = vld1q_u64(&sieve[w]);
        if ((vgetq_lane_u64(v,0)|vgetq_lane_u64(v,1))==0) continue;
        for (int k = 0; k < 2; ++k) {
            uint64_t word = sieve[w + k];
            while (word) {
                int bit = __builtin_ctzll(word);
                uint64_t idx = (w + k)*64 + bit;
                out_vec.push_back(static_cast<uint32_t>(low + 2 * idx));
                word &= word - 1;
            }
        }
    }
#endif
    for (; w < words; ++w) {
        uint64_t word = sieve[w];
        while (word) {
            int bit = __builtin_ctzll(word);
            uint64_t idx = w*64 + bit;
            out_vec.push_back(static_cast<uint32_t>(low + 2 * idx));
            word &= word - 1;
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " {slice_number}\n";
        return 1;
    }

    // Parse slice
    uint64_t slice = 0;
    try { slice = std::stoull(argv[1]); }
    catch (...) {
        std::cerr << "Invalid slice number\n";
        return 1;
    }

    // Build AA/BB/CC and DD paths
    unsigned b0 = slice>>24&0xFF,
             b1 = slice>>16&0xFF,
             b2 = slice>> 8&0xFF,
             b3 = slice    &0xFF;

    std::ostringstream dir_ss;
    dir_ss << std::uppercase<<std::hex<<std::setw(2)<<std::setfill('0')
           << b0<<"/"<<std::setw(2)<<b1<<"/"<<std::setw(2)<<b2;
    std::string dir = dir_ss.str();
    std::filesystem::create_directories(dir);

    std::ostringstream f_ss;
    f_ss << dir<<"/"
         << std::uppercase<<std::hex<<std::setw(2)<<std::setfill('0')
         << b3<<".dat";
    std::string out_path = f_ss.str();

    // Load 32-bit primes
    std::ifstream in32("uiprimes32.dat", std::ios::binary);
    if (!in32) {
        std::cerr << "Cannot open uiprimes32.dat\n";
        return 1;
    }
    in32.seekg(0, std::ios::end);
    size_t count32 = in32.tellg()/sizeof(uint32_t);
    in32.seekg(0);
    std::vector<uint32_t> small_primes(count32);
    in32.read(reinterpret_cast<char*>(small_primes.data()),
              count32*sizeof(uint32_t));
    in32.close();

    // Open slice output
    std::ofstream out(out_path, std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open " << out_path << "\n";
        return 1;
    }
    if (slice == 0) {
        uint32_t two = 2;
        out.write(reinterpret_cast<char*>(&two), sizeof(two));
    }

    // 64-bit interval
    uint64_t low64   = slice << 32;
    uint64_t high64  = ((slice+1)<<32) - 1;
    uint64_t low_odd = low64 | 1ULL;
    uint64_t total_odds = (high64 - low_odd)/2 + 1;

    // Segment size: 32 MiB → bits
    constexpr uint64_t SEG_BYTES = 32ULL*1024*1024;
    constexpr uint64_t SEG_BITS  = SEG_BYTES*8;
    uint64_t num_segs = (total_odds + SEG_BITS -1)/SEG_BITS;

    std::vector<uint64_t> sieve;
    std::vector<uint32_t> segment_primes;
    segment_primes.reserve(SEG_BITS/16);

    for (uint64_t seg = 0; seg < num_segs; ++seg) {
        uint64_t bit0  = seg * SEG_BITS;
        uint64_t bits  = std::min(SEG_BITS, total_odds - bit0);
        uint64_t start = low_odd + 2*bit0;
        uint64_t words = (bits + 63)/64;

        // Initialize all bits = “prime”
        sieve.assign(words, ~0ULL);
        if (bits & 63)
            sieve[words-1] = ~0ULL >> (64 - (bits & 63));

        // Clear “1” in slice=0, seg=0
        if (slice==0 && seg==0) sieve[0] &= ~1ULL;

        // Mark composites
        for (uint32_t prime : small_primes) {
            if (prime < 3) continue;
            uint64_t prime_sq = static_cast<uint64_t>(prime) * prime;
            if (prime_sq > high64) break;

            // Overflow-safe first multiple ≥ start
            __uint128_t rem  = static_cast<__uint128_t>(start) % prime;
            __uint128_t m128 = rem ? static_cast<__uint128_t>(start) + (prime - rem)
                                   : static_cast<__uint128_t>(start);
            if (m128 < prime_sq) m128 = prime_sq;
            if ((m128 & 1) == 0) m128 += prime;
            uint64_t bi = static_cast<uint64_t>((m128 - start) >> 1);

            for (uint64_t b = bi; b < bits; b += prime) {
                sieve[b>>6] &= ~(1ULL << (b & 63));
            }
        }

        // Collect & write
        segment_primes.clear();
        collect_primes_simd(sieve.data(), words, start, bits, segment_primes);
        out.write(reinterpret_cast<char*>(segment_primes.data()),
                  segment_primes.size() * sizeof(uint32_t));
    }

    out.close();
    std::cout << "Slice 0x" << std::hex << slice
              << " → " << out_path << "\n";
    return 0;
}
