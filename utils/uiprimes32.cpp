#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm> // For std::min, std::max
#include <string>    // For std::atoi
#include <cstdlib>   // For std::atoi
#include <thread>    // For std::thread::hardware_concurrency, std::thread
#include <future>    // For std::async, std::future
#include <vector>    // Already included, but for std::vector of futures
#include <functional> // For std::ref, std::cref if needed with std::bind or threads

// Intrinsics headers
#if defined(_MSC_VER)
#include <intrin.h>       // For MSVC intrinsics like _BitScanForward64, _mm_prefetch
#elif defined(__GNUC__) || defined(__clang__)
// GCC/Clang specific headers are usually included per function/architecture
#if defined(__x86_64__)
#include <immintrin.h>    // x86 AVX2
#elif defined(__aarch64__)
#include <arm_neon.h>     // ARM NEON
#endif
#else
// Potentially other compilers or no specific intrinsics
#endif

// --- Portability wrappers for builtins ---
// Count Trailing Zeros (ctzll)
#if defined(_MSC_VER)
inline int count_trailing_zeros_u64(unsigned long long val) {
    unsigned long index;
    if (_BitScanForward64(&index, val)) {
        return static_cast<int>(index);
    }
    return 64; // Standard behavior: undefined if val is 0; 64 is a common return
}
#elif defined(__GNUC__) || defined(__clang__)
#define count_trailing_zeros_u64 __builtin_ctzll
#else
// Fallback for other compilers (less efficient)
inline int count_trailing_zeros_u64(uint64_t n) {
    if (n == 0) return 64;
    int count = 0;
    while ((n & 1) == 0) {
        n >>= 1;
        count++;
        if (count >= 64) break; // Safety for n=0 if not caught above
    }
    return count;
}
#endif

// Prefetch
#if defined(_MSC_VER)
#define PREFETCH_WRITE(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_NTA)
#elif defined(__GNUC__) || defined(__clang__)
#define PREFETCH_WRITE(addr) __builtin_prefetch(addr, 1, 1) // rw=1 (write), locality=1 (low, NTA)
#else
#define PREFETCH_WRITE(addr) ((void)0) // No-op prefetch for other compilers
#endif
// --- End Portability wrappers ---


// SIMD-accelerated collection of primes from the bit-sieve
#if defined(__GNUC__) || defined(__clang__) // GCC/Clang attribute
#if defined(__x86_64__)
__attribute__((target("avx2")))
#endif
#endif
void collect_primes_simd(const uint64_t* sieve_data, uint64_t words,
                         uint64_t low_val_in_segment,
                         std::vector<uint32_t>& outvec) {
    uint64_t w = 0;

    // x86 AVX2 path: process 4 words (256 bits) at a time
#if (defined(__GNUC__) || defined(__clang__)) && defined(__AVX2__)
    for (; w + 4 <= words; w += 4) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&sieve_data[w]));
        if (_mm256_testz_si256(v, v)) continue; // Skip if all bits in 256-bit vector are zero
        for (int k = 0; k < 4; ++k) {
            uint64_t word = sieve_data[w + k];
            while (word) {
                int bit = count_trailing_zeros_u64(word);
                uint64_t idx_bit = (w + k) * 64 + bit;
                outvec.push_back(static_cast<uint32_t>(low_val_in_segment + 2 * idx_bit));
                word &= word - 1; // Clear the least significant bit set
            }
        }
    }
// ARM NEON path: process 2 words (128 bits) at a time
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__aarch64__) && defined(__ARM_NEON)
    for (; w + 2 <= words; w += 2) {
        uint64x2_t v = vld1q_u64(&sieve_data[w]);
        // if both lanes zero, skip
        if ((vgetq_lane_u64(v, 0) | vgetq_lane_u64(v, 1)) == 0) continue;
        for (int k = 0; k < 2; ++k) {
            uint64_t word = sieve_data[w + k];
            while (word) {
                int bit = count_trailing_zeros_u64(word);
                uint64_t idx_bit = (w + k) * 64 + bit;
                outvec.push_back(static_cast<uint32_t>(low_val_in_segment + 2 * idx_bit));
                word &= word - 1; // Clear the least significant bit set
            }
        }
    }
#elif defined(_MSC_VER) && defined(_M_X64) // MSVC AVX2 path (ensure /arch:AVX2 is set)
    // A proper runtime check for AVX2 support might involve __cpuidex (MSVC)
    // For simplicity, this path relies on compiler flags enabling AVX2.
    if (true) { // Placeholder for potential runtime AVX2 check
        for (; w + 4 <= words; w += 4) {
             __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&sieve_data[w]));
            if (_mm256_testz_si256(v, v)) continue;
            for (int k = 0; k < 4; ++k) {
                uint64_t word = sieve_data[w + k];
                while (word) {
                    int bit = count_trailing_zeros_u64(word);
                    uint64_t idx_bit = (w + k) * 64 + bit;
                    outvec.push_back(static_cast<uint32_t>(low_val_in_segment + 2 * idx_bit));
                    word &= word - 1;
                }
            }
        }
    }
#endif
    // Fallback: process remaining words one at a time (or if no SIMD path taken)
    for (; w < words; ++w) {
        uint64_t word = sieve_data[w];
        while (word) {
            int bit = count_trailing_zeros_u64(word);
            uint64_t idx_bit = w * 64 + bit;
            outvec.push_back(static_cast<uint32_t>(low_val_in_segment + 2 * idx_bit));
            word &= word - 1; // Clear the least significant bit set
        }
    }
}

// Function to process a single segment, designed to be run in a separate thread
std::vector<uint32_t> process_segment_task(
    uint32_t segment_idx_for_task,
    uint64_t total_odds_const,
    uint64_t segment_bits_capacity_const,
    const std::vector<uint32_t>& base_primes_const_ref) {

    // Determine the range of global bit indices for the current segment
    const uint64_t current_segment_low_bit_global_idx  = static_cast<uint64_t>(segment_idx_for_task) * segment_bits_capacity_const;
    const uint64_t current_segment_high_bit_global_idx = std::min(current_segment_low_bit_global_idx + segment_bits_capacity_const - 1, total_odds_const - 1);

    // Smallest odd number represented by the current segment
    const uint64_t low_val_in_segment = 3 + 2 * current_segment_low_bit_global_idx;
    // Largest odd number represented by the current segment
    const uint64_t high_val_in_segment = 3 + 2 * current_segment_high_bit_global_idx;

    // Number of actual bits used in this segment
    const uint64_t current_segment_active_bits = current_segment_high_bit_global_idx - current_segment_low_bit_global_idx + 1;
    // Number of 64-bit words needed for these active bits
    const uint64_t words_in_segment = (current_segment_active_bits + 63) / 64;

    std::vector<uint64_t> sieve_local(words_in_segment);
    std::vector<uint32_t> segment_primes_local;

    // Initialize bit-sieve for the current segment
    sieve_local.assign(words_in_segment, ~0ULL);
    if (current_segment_active_bits % 64 != 0) {
        sieve_local[words_in_segment - 1] &= (1ULL << (current_segment_active_bits % 64)) - 1;
    }

    // Mark composites in the current segment
    for (uint32_t p : base_primes_const_ref) {
        if (p == 2) continue;
        uint64_t p_squared = static_cast<uint64_t>(p) * p;
        if (p_squared > high_val_in_segment) break;

        uint64_t start_multiple = ((low_val_in_segment + p - 1) / p) * p;
        if (start_multiple < p_squared) start_multiple = p_squared;
        if ((start_multiple & 1) == 0) start_multiple += p;
        if (start_multiple > high_val_in_segment) continue;

        uint64_t bit_index = (start_multiple - low_val_in_segment) / 2;
        PREFETCH_WRITE(&sieve_local[bit_index >> 6]);
        for (uint64_t b = bit_index; b < current_segment_active_bits; b += p) {
            sieve_local[b >> 6] &= ~(1ULL << (b & 63));
        }
    }

    // Collect primes
    double density_approx_low = std::log(low_val_in_segment > 1 ? static_cast<double>(low_val_in_segment) : 2.0);
    if (density_approx_low < 0.1) density_approx_low = 0.1; // Prevent too small divisor, log(2) ~ 0.693
    segment_primes_local.reserve(static_cast<size_t>(static_cast<double>(current_segment_active_bits) / density_approx_low));

    collect_primes_simd(sieve_local.data(), words_in_segment, low_val_in_segment, segment_primes_local);
    return segment_primes_local;
}


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false); // Potentially speed up C++ streams
    std::cin.tie(NULL); // If using cin, but not critical here

    constexpr uint64_t MAX_N = 0xFFFFFFFFULL;
    int mb = 32; // Default segment size
    if (argc > 1) {
        mb = std::max(1, std::min(1024, std::atoi(argv[1])));
    }
    const uint64_t SEGMENT_BYTES = static_cast<uint64_t>(mb) * 1024 * 1024;
    const uint64_t segment_bits_capacity = SEGMENT_BYTES * 8;

    const auto sqrtN = static_cast<uint32_t>(std::sqrt(static_cast<double>(MAX_N)));
    std::vector<bool> is_prime_small(sqrtN + 1, true);
    is_prime_small[0] = is_prime_small[1] = false;
    std::vector<uint32_t> base_primes; // Primes up to sqrtN
    double reserve_val_log = std::log(sqrtN > 1 ? static_cast<double>(sqrtN) : 2.0);
    if (reserve_val_log < 0.1) reserve_val_log = 0.1;
    base_primes.reserve(static_cast<size_t>(static_cast<double>(sqrtN) / reserve_val_log));

    for (uint32_t i = 2; i <= sqrtN; ++i) {
        if (is_prime_small[i]) {
            base_primes.push_back(i);
            for (uint64_t j = static_cast<uint64_t>(i) * i; j <= sqrtN; j += i) {
                is_prime_small[j] = false;
            }
        }
    }

    const uint64_t effective_max_n_for_odds = (MAX_N % 2 == 0) ? MAX_N - 1 : MAX_N;
    const uint64_t total_odds = (effective_max_n_for_odds - 3) / 2 + 1;
    const auto num_segments = static_cast<uint32_t>((total_odds + segment_bits_capacity - 1) / segment_bits_capacity);

    std::ofstream out("uiprimes32.dat", std::ios::binary);
    if (!out) {
        std::cerr << "Error: Cannot open output file 'uiprimes32.dat'.\n";
        return 1;
    }
    uint32_t two = 2;
    out.write(reinterpret_cast<const char*>(&two), sizeof(two));

    // Determine number of threads to use
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1; // Fallback if detection fails
    // You might want to cap num_threads, e.g., std::min(num_threads, 8U) if too many threads cause overhead
    std::cout << "Using " << num_threads << " hardware threads.\n";
    std::cout << "Segment size: " << mb << "MB.\n";


    for (uint32_t batch_start_segment_idx = 0; batch_start_segment_idx < num_segments; batch_start_segment_idx += num_threads) {
        std::vector<std::future<std::vector<uint32_t>>> futures;
        futures.reserve(num_threads);

        // Launch a batch of segment processing tasks
        for (unsigned int i = 0; i < num_threads; ++i) {
            uint32_t current_segment_to_process = batch_start_segment_idx + i;
            if (current_segment_to_process >= num_segments) break; // Don't launch more tasks than segments

            futures.push_back(std::async(std::launch::async,
                                         process_segment_task,
                                         current_segment_to_process,
                                         total_odds,
                                         segment_bits_capacity,
                                         std::cref(base_primes))); // Pass base_primes by const reference
        }

        // Collect results from the batch and write to file in order
        for (size_t i = 0; i < futures.size(); ++i) {
            uint32_t actual_segment_idx_processed = batch_start_segment_idx + static_cast<uint32_t>(i);
            std::vector<uint32_t> primes_from_segment = futures[i].get(); // Waits for the task to complete

            if (!primes_from_segment.empty()) {
                out.write(reinterpret_cast<const char*>(primes_from_segment.data()),
                          static_cast<std::streamsize>(primes_from_segment.size() * sizeof(uint32_t)));
            }

            // Progress update for each completed segment in the batch
            const uint64_t current_segment_low_bit_global_idx  = static_cast<uint64_t>(actual_segment_idx_processed) * segment_bits_capacity;
            const uint64_t current_segment_high_bit_global_idx = std::min(current_segment_low_bit_global_idx + segment_bits_capacity - 1, total_odds - 1);
            const uint64_t high_val_in_segment_for_progress = 3 + 2 * current_segment_high_bit_global_idx;

            // More granular progress update
            if ((actual_segment_idx_processed +1) % 1 == 0 || actual_segment_idx_processed == num_segments -1 ) { // Update frequently or for last
                 std::cout << "Processed and wrote segment " << actual_segment_idx_processed + 1 << "/" << num_segments
                           << " (up to " << high_val_in_segment_for_progress << ")\n";
            }
        }
    }

    out.close();
    std::cout << "Prime generation complete. Primes stored in 'uiprimes32.dat'.\n";
    std::cout << "Configuration: Multi-threaded (" << num_threads << " threads), SIMD enabled (platform specific), "
              << mb << "MB segments.\n";
    return 0;
}

