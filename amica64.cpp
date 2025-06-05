#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <numeric>   // for std::gcd
#include <thread>
// #include <mutex> // No longer needed for pair output synchronization
#include <algorithm> // for std::min, std::sort
#include <sstream>   // for ostringstream

//------------------------------------------------------------------------------
// Global list of 32-bit primes
static std::vector<uint32_t> primes;
// static std::mutex io_mutex; // Removed as output is now centralized

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

    // Comparator for sorting by the smaller number (n)
    bool operator<(const AmicablePairOutput& other) const {
        if (n != other.n) {
            return n < other.n;
        }
        // If n is the same, sort by s (shouldn't happen for distinct amicable pairs where n < s)
        return s < other.s;
    }
};


// Load primes from a binary file
void load_primes(const std::string& filename) {
    std::ifstream in{filename, std::ios::binary};
    if (!in) {
        std::cerr << "Error: cannot open primes file '"
                  << filename << "'\n";
        std::exit(1);
    }
    uint32_t p_val;
    while (in.read(reinterpret_cast<char*>(&p_val), sizeof(p_val))) {
        primes.push_back(p_val);
    }
    std::cout << "Loaded " << primes.size()
              << " primes from " << filename << "\n";
}

// Calculates sum of proper divisors of n_val (sigma(n) - n)
uint64_t sum_proper_divisors(uint64_t n_val) {
    if (n_val < 2) return 0; // Proper divisors of 0 or 1 sum to 0

    uint64_t original_n = n_val;
    uint64_t sum_all_divs = 1;
    uint64_t temp_n = n_val;

    for (uint32_t p_u32 : primes) {
        uint64_t p = p_u32;
        if (p * p > temp_n) break;

        if (temp_n % p == 0) {
            uint64_t term_sum = 1;
            uint64_t p_power = 1;
            while (temp_n % p == 0) {
                temp_n /= p;
                p_power *= p;
                term_sum += p_power;
            }
            sum_all_divs *= term_sum;
        }
    }

    if (temp_n > 1) { // Remaining temp_n is a prime factor
        sum_all_divs *= (1 + temp_n);
    }
    return sum_all_divs - original_n;
}

// Factor x into (prime, exponent) list. Factors are sorted by prime.
std::vector<std::pair<uint64_t, uint32_t>> factor(uint64_t n_val) {
    std::vector<std::pair<uint64_t, uint32_t>> factors_list;
    if (n_val < 2) return factors_list; // Factorization of 0 or 1 is empty

    uint64_t temp_n = n_val;

    for (uint32_t p_u32 : primes) {
        uint64_t p = p_u32;
        if (p * p > temp_n) break;

        if (temp_n % p == 0) {
            uint32_t exponent = 0;
            while (temp_n % p == 0) {
                temp_n /= p;
                exponent++;
            }
            factors_list.emplace_back(p, exponent);
        }
    }

    if (temp_n > 1) { // Remaining temp_n is a prime factor
        factors_list.emplace_back(temp_n, 1);
    }
    return factors_list;
}


// Helper to derive prime factors of gcd(A, B) given sorted factors of A and B
std::vector<std::pair<uint64_t, uint32_t>> derive_gcd_factors(
    const std::vector<std::pair<uint64_t, uint32_t>>& factors1,
    const std::vector<std::pair<uint64_t, uint32_t>>& factors2) {

    std::vector<std::pair<uint64_t, uint32_t>> gcd_factors;
    auto it1 = factors1.begin();
    auto it2 = factors2.begin();

    while (it1 != factors1.end() && it2 != factors2.end()) {
        if (it1->first < it2->first) {
            ++it1;
        } else if (it2->first < it1->first) {
            ++it2;
        } else { // Same prime
            gcd_factors.emplace_back(it1->first, std::min(it1->second, it2->second));
            ++it1;
            ++it2;
        }
    }
    return gcd_factors;
}

// Helper to derive prime factors of num/divisor given their sorted factors
std::vector<std::pair<uint64_t, uint32_t>> derive_quotient_factors(
    const std::vector<std::pair<uint64_t, uint32_t>>& num_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& divisor_factors) {

    std::vector<std::pair<uint64_t, uint32_t>> quotient_factors;
    auto it_num = num_factors.begin();
    auto it_div = divisor_factors.begin();

    while (it_num != num_factors.end()) {
        if (it_div == divisor_factors.end() || it_num->first < it_div->first) {
            quotient_factors.push_back(*it_num);
            ++it_num;
        } else if (it_num->first > it_div->first) {
            ++it_div;
        } else {
            if (it_num->second > it_div->second) {
                quotient_factors.emplace_back(it_num->first, it_num->second - it_div->second);
            }
            ++it_num;
            ++it_div;
        }
    }
    return quotient_factors;
}

// Classifies an amicable pair (n, s)
ClassificationResult classify_amicable_pair(
    uint64_t n_val, uint64_t s_val,
    const std::vector<std::pair<uint64_t, uint32_t>>& n_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& s_factors) {

    ClassificationResult result;
    uint64_t g = std::gcd(n_val, s_val);

    auto g_factors = derive_gcd_factors(n_factors, s_factors);
    auto M_factors = derive_quotient_factors(n_factors, g_factors);
    auto N_factors = derive_quotient_factors(s_factors, g_factors);

    uint64_t M_val = n_val / g;
    uint64_t N_val = s_val / g;

    bool sqfM = true;
    for (const auto& factor_pair : M_factors) {
        if (factor_pair.second > 1) {
            sqfM = false;
            break;
        }
    }

    bool sqfN = true;
    for (const auto& factor_pair : N_factors) {
        if (factor_pair.second > 1) {
            sqfN = false;
            break;
        }
    }

    bool cpM = (std::gcd(M_val, g) == 1);
    bool cpN = (std::gcd(N_val, g) == 1);

    bool regular = sqfM && sqfN && cpM && cpN;
    std::string prefix = regular ? "" : "X";
    result.type_str = prefix + std::to_string(M_factors.size()) + "," + std::to_string(N_factors.size());

    return result;
}


// Format factor list as "p1^e1*p2^e2*…" using ostringstream
std::string format_factors(const std::vector<std::pair<uint64_t,uint32_t>>& factors) {
    std::ostringstream oss;
    for (size_t i = 0; i < factors.size(); ++i) {
        auto [p, e] = factors[i];
        oss << p;
        if (e > 1) {
            oss << '^' << e;
        }
        if (i + 1 < factors.size()) {
            oss << '*';
        }
    }
    return oss.str();
}

int main(const int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <max_n (uint64)> [primes_file (optional)]\n";
        return 1;
    }
    uint64_t max_n_arg = 0;
    try {
        max_n_arg = std::stoull(argv[1]);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: max_n value '" << argv[1] << "' is out of range for uint64_t.\n";
        return 1;
    } catch (const std::invalid_argument& ia) {
        std::cerr << "Error: max_n value '" << argv[1] << "' is not a valid unsigned integer.\n";
        return 1;
    }

    if (max_n_arg < 220) {
        std::cout << "max_n (" << max_n_arg << ") is too small to find amicable pairs. Minimum is 220.\n";
        return 0;
    }

    const std::string primes_file_path = argc == 3 ? argv[2] : "uiprimes32.dat";
    load_primes(primes_file_path);
    if (primes.empty()) {
        std::cerr << "Error: No primes loaded. Ensure '" << primes_file_path
                  << "' is valid, accessible, and contains primes.\n";
        return 1;
    }

    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) {
        std::cerr << "Error: cannot open output file amicable_pairs.txt\n";
        return 1;
    }

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 2;

    const uint64_t search_limit_val = max_n_arg;
    const uint64_t total_numbers_to_check = (search_limit_val >= 2) ? (search_limit_val - 1) : 0;


    if (total_numbers_to_check == 0) {
        std::cout << "No numbers in the specified range (2 to " << search_limit_val << ") to check.\n";
        if (out_file_stream.is_open()) out_file_stream.close();
        return 0;
    }

    uint64_t chunk_size = total_numbers_to_check / num_threads;
    if (chunk_size == 0 && total_numbers_to_check > 0) chunk_size = 1;
    if (total_numbers_to_check < num_threads) num_threads = static_cast<unsigned int>(total_numbers_to_check);


    std::vector<std::thread> threads;
    threads.reserve(num_threads);

    std::vector<std::vector<AmicablePairOutput>> per_thread_results(num_threads);


    std::cout << "Searching for amicable pairs up to " << search_limit_val
              << " using " << num_threads << " threads.\n";

    for (unsigned int t_idx = 0; t_idx < num_threads; ++t_idx) {
        uint64_t start_n_val = 2 + t_idx * chunk_size;
        uint64_t end_n_val = (t_idx == num_threads - 1) ? search_limit_val : (start_n_val + chunk_size - 1);

        if (start_n_val > search_limit_val) continue;

        // Pass a reference to this thread's result vector
        threads.emplace_back([start_n_val, end_n_val, search_limit_val, &thread_specific_results = per_thread_results[t_idx]] {
            for (uint64_t n_current = start_n_val; n_current <= end_n_val; ++n_current) {
                if (n_current < 220) { // Smallest n in amicable pair (n,s) with n < s is 220
                    continue;
                }

                uint64_t s_candidate = sum_proper_divisors(n_current);

                if (s_candidate > n_current && s_candidate <= search_limit_val) {
                    if (sum_proper_divisors(s_candidate) == n_current) {
                        // Amicable pair found, factor and classify
                        auto n_factors = factor(n_current);
                        auto s_factors = factor(s_candidate);

                        ClassificationResult classification = classify_amicable_pair(
                            n_current, s_candidate,
                            n_factors, s_factors
                        );

                        std::string fn_str = format_factors(n_factors);
                        std::string fs_str = format_factors(s_factors);

                        AmicablePairOutput pair_data;
                        pair_data.n = n_current;
                        pair_data.s = s_candidate;
                        pair_data.classification_str = classification.type_str;
                        pair_data.n_factors_str = fn_str;
                        pair_data.s_factors_str = fs_str;
                        thread_specific_results.push_back(pair_data);
                    }
                }
            }
        });
    }

    for (auto &th : threads) {
        if (th.joinable()) {
            th.join();
        }
    }

    // Consolidate results from all threads
    std::vector<AmicablePairOutput> all_found_pairs;
    for (const auto& thread_results : per_thread_results) {
        all_found_pairs.insert(all_found_pairs.end(), thread_results.begin(), thread_results.end());
    }

    // Sort all found pairs
    std::sort(all_found_pairs.begin(), all_found_pairs.end());

    // Print sorted results
    for (const auto& pair_data : all_found_pairs) {
        out_file_stream << pair_data.classification_str << '\n'
                        << pair_data.n  << '=' << pair_data.n_factors_str << '\n'
                        << pair_data.s  << '=' << pair_data.s_factors_str << "\n\n";
        std::cout << pair_data.classification_str << '\n'
                  << pair_data.n  << '=' << pair_data.n_factors_str << '\n'
                  << pair_data.s  << '=' << pair_data.s_factors_str << "\n\n";
    }

    std::cout << "Done. Results in amicable_pairs.txt\n";
    if (out_file_stream.is_open()) {
        out_file_stream.close();
    }
    return 0;
}

