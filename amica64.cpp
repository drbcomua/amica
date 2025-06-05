#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <numeric>   // for std::gcd
#include <thread>
#include <algorithm> // for std::min, std::sort, std::max
#include <sstream>   // for ostringstream
#include <cmath>     // for std::sqrt, std::log
#include <stdexcept> // for std::stoi/stoull exceptions
#include <cstdlib>   // For std::atoi in prime gen (can be replaced with std::stoi)
#include <future>    // For std::async, std::future
#include <functional> // For std::ref, std::cref

// SIMD Intrinsics Headers
#if defined(_MSC_VER)
#include <intrin.h>       // For MSVC intrinsics like _BitScanForward64, _mm_prefetch
#elif defined(__GNUC__) || defined(__clang__)
// GCC/Clang specific headers are usually included per function/architecture
#if defined(__x86_64__) || defined(_M_X64) // Added _M_X64 for MSVC x64 compatibility
#include <immintrin.h>    // x86 AVX2
#elif defined(__aarch64__)
#include <arm_neon.h>     // ARM NEON
#endif
#endif


//------------------------------------------------------------------------------
// Portability wrappers for builtins
//------------------------------------------------------------------------------
#if defined(_MSC_VER)
inline int count_trailing_zeros_u64(unsigned long long val) {
    unsigned long index;
    if (_BitScanForward64(&index, val)) {
        return static_cast<int>(index);
    }
    return 64;
}
#define PREFETCH_WRITE(addr) _mm_prefetch(reinterpret_cast<const char*>(addr), _MM_HINT_NTA)
#elif defined(__GNUC__) || defined(__clang__)
#define count_trailing_zeros_u64 __builtin_ctzll
#define PREFETCH_WRITE(addr) __builtin_prefetch(addr, 1, 1) // rw=1 (write), locality=1 (low, NTA)
#else
// Fallback for other compilers (less efficient)
inline int count_trailing_zeros_u64(uint64_t n) {
    if (n == 0) return 64;
    int count = 0;
    while ((n & 1) == 0) {
        n >>= 1;
        count++;
        if (count >= 64) break;
    }
    return count;
}
#define PREFETCH_WRITE(addr) ((void)0) // No-op prefetch
#endif

//------------------------------------------------------------------------------
// Global list of 32-bit primes
static std::vector<uint32_t> primes;

// Struct to hold the result of amicable pair classification
struct ClassificationResult {
    std::string type_str;
};

// Struct to hold all data for an amicable pair for sorted output
struct AmicablePairOutput {
    uint64_t n;
    uint64_t s;
    std::string classification_str;
    std::string n_factors_str;
    std::string s_factors_str;

    bool operator<(const AmicablePairOutput& other) const {
        if (n != other.n) { return n < other.n; }
        return s < other.s;
    }
};

//------------------------------------------------------------------------------
// Prime Generation Sieve (Multi-Threaded)
//------------------------------------------------------------------------------

#if defined(__GNUC__) || defined(__clang__) // GCC/Clang attribute
#if defined(__x86_64__)
__attribute__((target("avx2")))
#endif
#endif
void collect_primes_simd(const uint64_t* sieve_data, uint64_t words,
                         uint64_t low_val_in_segment,
                         std::vector<uint32_t>& outvec) {
    uint64_t w = 0;

#if (defined(__GNUC__) || defined(__clang__)) && defined(__AVX2__)
    for (; w + 4 <= words; w += 4) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&sieve_data[w]));
        if (_mm256_testz_si256(v, v)) continue;
        for (int k = 0; k < 4; ++k) {
            uint64_t word = sieve_data[w + k];
            while (word) {
                int bit = count_trailing_zeros_u64(word);
                uint64_t idx_bit = (w + k) * 64 + bit;
                uint64_t prime_candidate = low_val_in_segment + 2 * idx_bit;
                if (prime_candidate <= 0xFFFFFFFFULL) {
                     outvec.push_back(static_cast<uint32_t>(prime_candidate));
                }
                word &= word - 1;
            }
        }
    }
#elif (defined(__GNUC__) || defined(__clang__)) && defined(__aarch64__) && defined(__ARM_NEON)
    for (; w + 2 <= words; w += 2) {
        uint64x2_t v = vld1q_u64(&sieve_data[w]);
        if ((vgetq_lane_u64(v, 0) | vgetq_lane_u64(v, 1)) == 0) continue;
        for (int k = 0; k < 2; ++k) {
            uint64_t word = sieve_data[w + k];
            while (word) {
                int bit = count_trailing_zeros_u64(word);
                uint64_t idx_bit = (w + k) * 64 + bit;
                uint64_t prime_candidate = low_val_in_segment + 2 * idx_bit;
                 if (prime_candidate <= 0xFFFFFFFFULL) {
                    outvec.push_back(static_cast<uint32_t>(prime_candidate));
                }
                word &= word - 1;
            }
        }
    }
#elif defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86)) // MSVC AVX2 path (check __AVX2__ defined by compiler if /arch:AVX2)
    // A proper runtime check for AVX2 support might involve __cpuidex (MSVC)
    // For simplicity, this path relies on compiler flags enabling AVX2 if __AVX2__ is defined by compiler.
    #ifdef __AVX2__
    bool has_avx2 = true; // Assuming compiler defines __AVX2__ if available and enabled
    if (has_avx2) {
        for (; w + 4 <= words; w += 4) {
             __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&sieve_data[w]));
            if (_mm256_testz_si256(v, v)) continue;
            for (int k = 0; k < 4; ++k) {
                uint64_t word = sieve_data[w + k];
                while (word) {
                    int bit = count_trailing_zeros_u64(word);
                    uint64_t idx_bit = (w + k) * 64 + bit;
                    uint64_t prime_candidate = low_val_in_segment + 2 * idx_bit;
                    if (prime_candidate <= 0xFFFFFFFFULL) {
                        outvec.push_back(static_cast<uint32_t>(prime_candidate));
                    }
                    word &= word - 1;
                }
            }
        }
    }
    #endif // __AVX2__
#endif
    for (; w < words; ++w) {
        uint64_t word = sieve_data[w];
        while (word) {
            int bit = count_trailing_zeros_u64(word);
            uint64_t idx_bit = w * 64 + bit;
            uint64_t prime_candidate = low_val_in_segment + 2 * idx_bit;
            if (prime_candidate <= 0xFFFFFFFFULL) {
                outvec.push_back(static_cast<uint32_t>(prime_candidate));
            }
            word &= word - 1;
        }
    }
}

std::vector<uint32_t> process_segment_task(
    uint32_t segment_idx_for_task,
    uint64_t total_sieve_slots_count_const, // Renamed to avoid conflict
    uint64_t segment_bits_capacity_const,
    const std::vector<uint32_t>& base_primes_const_ref) {

    const uint64_t current_segment_low_bit_global_idx  = static_cast<uint64_t>(segment_idx_for_task) * segment_bits_capacity_const;
    if (current_segment_low_bit_global_idx >= total_sieve_slots_count_const) {
        return {}; // No work for this segment index
    }
    const uint64_t current_segment_high_bit_global_idx = std::min(current_segment_low_bit_global_idx + segment_bits_capacity_const - 1, total_sieve_slots_count_const - 1);

    const uint64_t low_val_in_segment = 3 + 2 * current_segment_low_bit_global_idx;
    const uint64_t high_val_in_segment = 3 + 2 * current_segment_high_bit_global_idx;

    const uint64_t current_segment_active_bits = current_segment_high_bit_global_idx - current_segment_low_bit_global_idx + 1;
    const uint64_t words_in_segment = (current_segment_active_bits + 63) / 64;

    std::vector<uint64_t> sieve_local(words_in_segment);
    std::vector<uint32_t> segment_primes_local;

    sieve_local.assign(words_in_segment, ~0ULL);
    if (current_segment_active_bits % 64 != 0) {
        sieve_local[words_in_segment - 1] &= (1ULL << (current_segment_active_bits % 64)) - 1;
    }

    for (uint32_t p : base_primes_const_ref) {
        if (p == 2) continue;
        uint64_t p_squared = static_cast<uint64_t>(p) * p;
        if (p_squared > high_val_in_segment) break;

        uint64_t start_multiple = ((low_val_in_segment + p - 1) / p) * p;
        if (start_multiple < p_squared) start_multiple = p_squared;
        if ((start_multiple & 1) == 0) start_multiple += p; // Ensure odd multiple
        if (start_multiple > high_val_in_segment) continue;

        uint64_t bit_index = (start_multiple - low_val_in_segment) / 2;
        PREFETCH_WRITE(&sieve_local[bit_index >> 6]);
        for (uint64_t b = bit_index; b < current_segment_active_bits; b += p) {
            sieve_local[b >> 6] &= ~(1ULL << (b & 63));
        }
    }

    double log_low_val = std::log(static_cast<double>(low_val_in_segment > 1 ? low_val_in_segment : 2.0));
    if (log_low_val < 0.1) log_low_val = 0.1; // Avoid division by zero or tiny number
    segment_primes_local.reserve(static_cast<size_t>(static_cast<double>(current_segment_active_bits) / log_low_val * 1.2) + 100); // 1.2 factor for buffer

    collect_primes_simd(sieve_local.data(), words_in_segment, low_val_in_segment, segment_primes_local);
    return segment_primes_local;
}

void generate_primes_multi_threaded(int mb_segment_size_param) {
    ::primes.clear();
    std::ios_base::sync_with_stdio(false); // Might not have a huge effect here but good practice
    std::cin.tie(NULL);

    constexpr uint64_t MAX_SIEVE_N = 0xFFFFFFFFULL;
    const uint64_t SEGMENT_BYTES = static_cast<uint64_t>(mb_segment_size_param) * 1024 * 1024;
    const uint64_t segment_bits_capacity = SEGMENT_BYTES * 8;

    const auto sqrt_max_sieve_n = static_cast<uint32_t>(std::sqrt(static_cast<double>(MAX_SIEVE_N)));
    std::vector<bool> is_small_prime(sqrt_max_sieve_n + 1, true);
    is_small_prime[0] = is_small_prime[1] = false;
    std::vector<uint32_t> base_primes;
    double log_sqrt_n = std::log(static_cast<double>(sqrt_max_sieve_n > 1 ? sqrt_max_sieve_n : 2.0));
    if (log_sqrt_n < 0.1) log_sqrt_n = 0.1;
    base_primes.reserve(static_cast<size_t>(static_cast<double>(sqrt_max_sieve_n) / log_sqrt_n * 1.2) + 100);


    for (uint32_t i = 2; i <= sqrt_max_sieve_n; ++i) {
        if (is_small_prime[i]) {
            base_primes.push_back(i);
            for (uint64_t j = static_cast<uint64_t>(i) * i; j <= sqrt_max_sieve_n; j += i) {
                is_small_prime[j] = false;
            }
        }
    }

    // total_sieve_slots_count: number of odd numbers from 3 to MAX_SIEVE_N
    // For MAX_SIEVE_N = 7 (3,5,7), count = 3. (7-1)/2 = 3.
    // For MAX_SIEVE_N = 8 (3,5,7), count = 3. (8-1)/2 = 3 (integer division).
    // Max index is count -1. Max value is 3 + 2*(count-1).
    const uint64_t total_sieve_slots_count = (MAX_SIEVE_N -1) / 2;

    const auto num_segments = static_cast<uint32_t>((total_sieve_slots_count + segment_bits_capacity - 1) / segment_bits_capacity);

    ::primes.push_back(2);
    // Estimate final prime count for reservation (approx pi(x) for x = MAX_SIEVE_N)
    double log_max_n = std::log(static_cast<double>(MAX_SIEVE_N > 1 ? MAX_SIEVE_N : 2.0));
    if (log_max_n < 0.1) log_max_n = 0.1;
    ::primes.reserve(static_cast<size_t>(static_cast<double>(MAX_SIEVE_N) / log_max_n * 1.1) + 100);


    unsigned int num_hw_threads = std::thread::hardware_concurrency();
    if (num_hw_threads == 0) num_hw_threads = 1;
    // Cap threads to avoid excessive overhead, e.g., max 16 or 32, or relative to num_segments
    num_hw_threads = std::min(num_hw_threads, static_cast<unsigned int>(32)); // Example cap
    num_hw_threads = std::min(num_hw_threads, num_segments > 0 ? num_segments : 1U); // Don't use more threads than segments


    std::cout << "Generating primes in memory up to " << MAX_SIEVE_N
              << " (using " << mb_segment_size_param << "MB segments, "
              << num_hw_threads << " threads)...\n";

    std::vector<std::vector<uint32_t>> results_from_all_segments(num_segments); // Store results per segment to ensure order

    for (uint32_t batch_start_segment_idx = 0; batch_start_segment_idx < num_segments; batch_start_segment_idx += num_hw_threads) {
        std::vector<std::future<std::vector<uint32_t>>> futures_in_batch;
        futures_in_batch.reserve(num_hw_threads);

        for (unsigned int i = 0; i < num_hw_threads; ++i) {
            uint32_t current_segment_to_process = batch_start_segment_idx + i;
            if (current_segment_to_process >= num_segments) break;

            futures_in_batch.push_back(std::async(std::launch::async,
                                         process_segment_task,
                                         current_segment_to_process,
                                         total_sieve_slots_count,
                                         segment_bits_capacity,
                                         std::cref(base_primes)));
        }

        for (size_t i = 0; i < futures_in_batch.size(); ++i) {
            uint32_t actual_segment_idx_processed = batch_start_segment_idx + static_cast<uint32_t>(i);
            results_from_all_segments[actual_segment_idx_processed] = futures_in_batch[i].get();

            if (((actual_segment_idx_processed + 1) % (num_segments / 100 + 1) == 0) || (actual_segment_idx_processed == num_segments - 1) ) {
                 std::cout << "\rProcessed segment " << actual_segment_idx_processed + 1 << "/" << num_segments << std::flush;
            }
        }
    }
    std::cout << "\nAll segments processed. Consolidating primes...\n";

    // Consolidate primes from results_from_all_segments into global ::primes
    for(const auto& segment_primes_vec : results_from_all_segments) {
        ::primes.insert(::primes.end(), segment_primes_vec.begin(), segment_primes_vec.end());
    }

    std::cout << "Prime generation complete. Found " << ::primes.size() << " primes.\n";
}

//------------------------------------------------------------------------------
// Amicable Pair Logic (largely unchanged)
//------------------------------------------------------------------------------

uint64_t sum_proper_divisors(uint64_t n_val) {
    if (n_val < 2) return 0;
    uint64_t original_n = n_val;
    uint64_t sum_all_divs = 1;
    uint64_t temp_n = n_val;
    for (uint32_t p_u32 : ::primes) {
        uint64_t p = p_u32;
        if (p == 0) continue;
        if (p * p > temp_n) break;
        if (temp_n % p == 0) {
            uint64_t term_sum = 1;
            uint64_t p_power = 1;
            while (temp_n % p == 0) { temp_n /= p; p_power *= p; term_sum += p_power; }
            sum_all_divs *= term_sum;
        }
    }
    if (temp_n > 1) { sum_all_divs *= (1 + temp_n); }
    return sum_all_divs - original_n;
}

std::vector<std::pair<uint64_t, uint32_t>> factor(uint64_t n_val) {
    std::vector<std::pair<uint64_t, uint32_t>> factors_list;
    if (n_val < 2) return factors_list;
    uint64_t temp_n = n_val;
    for (uint32_t p_u32 : ::primes) {
        uint64_t p = p_u32;
        if (p == 0) continue;
        if (p * p > temp_n) break;
        if (temp_n % p == 0) {
            uint32_t exponent = 0;
            while (temp_n % p == 0) { temp_n /= p; exponent++; }
            factors_list.emplace_back(p, exponent);
        }
    }
    if (temp_n > 1) { factors_list.emplace_back(temp_n, 1); }
    return factors_list;
}

std::vector<std::pair<uint64_t, uint32_t>> derive_gcd_factors(
    const std::vector<std::pair<uint64_t, uint32_t>>& factors1,
    const std::vector<std::pair<uint64_t, uint32_t>>& factors2) {
    std::vector<std::pair<uint64_t, uint32_t>> gcd_factors;
    auto it1 = factors1.begin(); auto it2 = factors2.begin();
    while (it1 != factors1.end() && it2 != factors2.end()) {
        if (it1->first < it2->first) { ++it1; }
        else if (it2->first < it1->first) { ++it2; }
        else { gcd_factors.emplace_back(it1->first, std::min(it1->second, it2->second)); ++it1; ++it2;}
    }
    return gcd_factors;
}

std::vector<std::pair<uint64_t, uint32_t>> derive_quotient_factors(
    const std::vector<std::pair<uint64_t, uint32_t>>& num_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& divisor_factors) {
    std::vector<std::pair<uint64_t, uint32_t>> quotient_factors;
    auto it_num = num_factors.begin(); auto it_div = divisor_factors.begin();
    while (it_num != num_factors.end()) {
        if (it_div == divisor_factors.end() || it_num->first < it_div->first) { quotient_factors.push_back(*it_num); ++it_num;}
        else if (it_num->first > it_div->first) { ++it_div; }
        else { if (it_num->second > it_div->second) { quotient_factors.emplace_back(it_num->first, it_num->second - it_div->second); } ++it_num; ++it_div; }
    }
    return quotient_factors;
}

ClassificationResult classify_amicable_pair(
    uint64_t n_val, uint64_t s_val,
    const std::vector<std::pair<uint64_t, uint32_t>>& n_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& s_factors) {
    ClassificationResult result;
    uint64_t g = std::gcd(n_val, s_val);
    auto g_factors = derive_gcd_factors(n_factors, s_factors);
    auto M_factors = derive_quotient_factors(n_factors, g_factors);
    auto N_factors = derive_quotient_factors(s_factors, g_factors);
    uint64_t M_val = n_val / g; uint64_t N_val = s_val / g;
    bool sqfM = true; for (const auto& f : M_factors) if (f.second > 1) { sqfM = false; break; }
    bool sqfN = true; for (const auto& f : N_factors) if (f.second > 1) { sqfN = false; break; }
    bool cpM = (std::gcd(M_val, g) == 1); bool cpN = (std::gcd(N_val, g) == 1);
    bool regular = sqfM && sqfN && cpM && cpN;
    result.type_str = (regular ? "" : "X") + std::to_string(M_factors.size()) + "," + std::to_string(N_factors.size());
    return result;
}

std::string format_factors(const std::vector<std::pair<uint64_t,uint32_t>>& factors) {
    std::ostringstream oss;
    for (size_t i = 0; i < factors.size(); ++i) {
        auto [p, e] = factors[i]; oss << p;
        if (e > 1) { oss << '^' << e; }
        if (i + 1 < factors.size()) { oss << '*'; }
    }
    return oss.str();
}

int main(const int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <max_n (uint64)> [segment_size_mb (optional, int, default 16)]\n";
        return 1;
    }
    uint64_t max_n_arg = 0;
    try {
        max_n_arg = std::stoull(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid max_n value '" << argv[1] << "'. " << e.what() << "\n";
        return 1;
    }

    if (max_n_arg < 220) {
        std::cout << "max_n (" << max_n_arg << ") is too small. Minimum is 220.\n";
        return 0;
    }

    int segment_size_mb_val = 16;
    if (argc == 3) {
        try {
            segment_size_mb_val = std::stoi(argv[2]);
            if (segment_size_mb_val < 1 || segment_size_mb_val > 2048) {
                std::cerr << "Warning: Segment size MB (" << segment_size_mb_val
                          << ") out of reasonable range [1, 2048]. Using default 16MB.\n";
                segment_size_mb_val = 16;
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Invalid segment size '" << argv[2]
                      << "'. Using default 16MB. Error: " << e.what() << "\n";
            segment_size_mb_val = 16;
        }
    }

    generate_primes_multi_threaded(segment_size_mb_val); // Changed function name
    if (::primes.empty() || ::primes.front() != 2) {
        std::cerr << "Error: Prime generation failed.\n";
        return 1;
    }

    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) {
        std::cerr << "Error: cannot open output file amicable_pairs.txt\n";
        return 1;
    }

    unsigned int num_threads_amicable = std::thread::hardware_concurrency(); // Different variable for amicable search threads
    if (num_threads_amicable == 0) num_threads_amicable = 2;
    num_threads_amicable = std::min(num_threads_amicable, static_cast<unsigned int>(64)); // Cap amicable search threads too

    const uint64_t search_limit_val = max_n_arg;
    const uint64_t first_n_to_check_loop = 2;
    const uint64_t total_numbers_for_chunking = (search_limit_val >= first_n_to_check_loop) ? (search_limit_val - first_n_to_check_loop + 1) : 0;

    if (total_numbers_for_chunking == 0) {
        std::cout << "No numbers in range to check.\n";
        if (out_file_stream.is_open()) out_file_stream.close();
        return 0;
    }

    uint64_t chunk_size = total_numbers_for_chunking / num_threads_amicable;
    if (chunk_size == 0 && total_numbers_for_chunking > 0) chunk_size = 1;
    if (total_numbers_for_chunking < num_threads_amicable) num_threads_amicable = static_cast<unsigned int>(total_numbers_for_chunking);
    if (num_threads_amicable == 0 && total_numbers_for_chunking > 0) num_threads_amicable = 1;

    std::vector<std::thread> threads_amicable_search; // Renamed vector
    threads_amicable_search.reserve(num_threads_amicable);
    std::vector<std::vector<AmicablePairOutput>> per_thread_results(num_threads_amicable);

    std::cout << "Searching for amicable pairs up to " << search_limit_val
              << " using " << num_threads_amicable << " threads.\n";

    for (unsigned int t_idx = 0; t_idx < num_threads_amicable; ++t_idx) {
        uint64_t start_n_val = first_n_to_check_loop + t_idx * chunk_size;
        uint64_t end_n_val = (t_idx == num_threads_amicable - 1) ? search_limit_val : (start_n_val + chunk_size - 1);
        if (start_n_val > search_limit_val) continue;

        threads_amicable_search.emplace_back([start_n_val, end_n_val, search_limit_val, &res_vec = per_thread_results[t_idx]] {
            for (uint64_t n = start_n_val; n <= end_n_val; ++n) {
                if (n < 220) continue;
                uint64_t s = sum_proper_divisors(n);
                if (s > n && s <= search_limit_val) {
                    if (sum_proper_divisors(s) == n) {
                        auto nf = factor(n); auto sf = factor(s);
                        auto cls = classify_amicable_pair(n, s, nf, sf);
                        res_vec.push_back({n, s, cls.type_str, format_factors(nf), format_factors(sf)});
                    }
                }
            }
        });
    }

    for (auto &th : threads_amicable_search) if (th.joinable()) th.join();

    std::vector<AmicablePairOutput> all_found_pairs;
    for (const auto& R : per_thread_results) all_found_pairs.insert(all_found_pairs.end(), R.begin(), R.end());
    std::sort(all_found_pairs.begin(), all_found_pairs.end());

    for (const auto& p : all_found_pairs) {
        out_file_stream << p.classification_str << '\n' << p.n  << '=' << p.n_factors_str << '\n' << p.s  << '=' << p.s_factors_str << "\n\n";
        std::cout       << p.classification_str << '\n' << p.n  << '=' << p.n_factors_str << '\n' << p.s  << '=' << p.s_factors_str << "\n\n";
    }

    std::cout << "Done. Results in amicable_pairs.txt\n";
    if (out_file_stream.is_open()) out_file_stream.close();
    return 0;
}

