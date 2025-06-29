#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <numeric>
#include <thread>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <future>
#include <functional>

//------------------------------------------------------------------------------
// Global Data & Structures
//------------------------------------------------------------------------------

static std::vector<uint32_t> primes;

struct AmicablePairOutput {
    uint64_t n;
    uint64_t s;
    std::string classification_str;
    std::string n_factors_str;
    std::string s_factors_str;

    bool operator<(const AmicablePairOutput& other) const {
        return n < other.n;
    }
};

//------------------------------------------------------------------------------
// Prime Generation (Optimized and Scalable)
//------------------------------------------------------------------------------
void generate_primes(uint64_t limit) {
    ::primes.clear();
    std::cout << "Generating primes up to " << limit << " for factorization..." << std::endl;
    if (limit < 2) return;

    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (uint64_t p = 2; p * p <= limit; ++p) {
        if (is_prime[p]) {
            for (uint64_t i = p * p; i <= limit; i += p)
                is_prime[i] = false;
        }
    }
    ::primes.push_back(2);
    for (uint64_t p = 3; p <= limit; p += 2) {
        if (is_prime[p]) {
            ::primes.push_back(static_cast<uint32_t>(p));
        }
    }
    std::cout << "Prime generation complete. Found " << ::primes.size() << " primes." << std::endl;
}

//------------------------------------------------------------------------------
// Amicable Pair Logic
//------------------------------------------------------------------------------

uint64_t sum_proper_divisors(uint64_t n_val) {
    if (n_val < 2) return 0;
    uint64_t original_n = n_val;
    uint64_t sum_all_divs = 1;
    uint64_t temp_n = n_val;
    for (uint32_t p_u32 : ::primes) {
        uint64_t p = p_u32;
        if (p * p > temp_n) break;
        if (temp_n % p == 0) {
            uint64_t term_sum = 1;
            uint64_t p_power = 1;
            do {
                temp_n /= p;
                p_power *= p;
                term_sum += p_power;
            } while (temp_n % p == 0);
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

std::string classify_amicable_pair(
    uint64_t n_val, uint64_t s_val,
    const std::vector<std::pair<uint64_t, uint32_t>>& n_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& s_factors) {
    uint64_t g = std::gcd(n_val, s_val);
    auto g_factors = derive_gcd_factors(n_factors, s_factors);
    auto M_factors = derive_quotient_factors(n_factors, g_factors);
    auto N_factors = derive_quotient_factors(s_factors, g_factors);
    uint64_t M_val = n_val / g; uint64_t N_val = s_val / g;
    bool sqfM = true; for (const auto& f : M_factors) if (f.second > 1) { sqfM = false; break; }
    bool sqfN = true; for (const auto& f : N_factors) if (f.second > 1) { sqfN = false; break; }
    bool cpM = (std::gcd(M_val, g) == 1); bool cpN = (std::gcd(N_val, g) == 1);
    bool regular = sqfM && sqfN && cpM && cpN;
    return (regular ? "" : "X") + std::to_string(M_factors.size()) + "," + std::to_string(N_factors.size());
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

//------------------------------------------------------------------------------
// Main Program Logic
//------------------------------------------------------------------------------
int main(const int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <max_n (uint64)>\n";
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

    // Generate primes up to a sensible limit. s(n) can be larger than n, so we give
    // a buffer. sqrt(max_n * 2) is a reasonable heuristic.
    double safe_upper_bound = static_cast<double>(max_n_arg) * 2.0;
    uint64_t prime_limit = static_cast<uint64_t>(std::sqrt(safe_upper_bound)) + 1;
    generate_primes(prime_limit);

    if (::primes.empty()) {
        std::cerr << "Error: Prime generation failed.\n";
        return 1;
    }

    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) {
        std::cerr << "Error: cannot open output file amicable_pairs.txt\n";
        return 1;
    }

    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 2;
    std::cout << "Searching for amicable pairs up to " << max_n_arg
              << " using " << num_threads << " threads.\n";

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    std::vector<std::vector<AmicablePairOutput>> per_thread_results(num_threads);

    uint64_t chunk_size = max_n_arg / num_threads;
    if (chunk_size < 1) chunk_size = 1;

    for (unsigned int i = 0; i < num_threads; ++i) {
        uint64_t start_n = i * chunk_size + 1;
        uint64_t end_n = (i == num_threads - 1) ? max_n_arg : (start_n + chunk_size - 1);
        if (start_n > max_n_arg) continue;

        threads.emplace_back([start_n, end_n, max_n_arg, &res_vec = per_thread_results[i]] {
            for (uint64_t n = start_n; n <= end_n; ++n) {
                if (n < 220) continue;
                uint64_t s = sum_proper_divisors(n);
                if (s > n && s <= max_n_arg) {
                    if (sum_proper_divisors(s) == n) {
                        auto nf = factor(n);
                        auto sf = factor(s);
                        auto cls = classify_amicable_pair(n, s, nf, sf);
                        res_vec.push_back({n, s, cls, format_factors(nf), format_factors(sf)});
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

    std::vector<AmicablePairOutput> all_found_pairs;
    for (const auto& result_chunk : per_thread_results) {
        all_found_pairs.insert(all_found_pairs.end(), result_chunk.begin(), result_chunk.end());
    }
    std::sort(all_found_pairs.begin(), all_found_pairs.end());

    std::cout << "\n--- Found " << all_found_pairs.size() << " Amicable Pairs ---\n\n";
    for (const auto& p : all_found_pairs) {
        std::ostringstream oss;
        oss << p.classification_str << '\n' << p.n  << '=' << p.n_factors_str << '\n' << p.s  << '=' << p.s_factors_str << "\n\n";
        std::cout << oss.str();
        out_file_stream << oss.str();
    }

    std::cout << "Done. Found a total of " << all_found_pairs.size() << " pairs. Results also saved in amicable_pairs.txt\n";
    if (out_file_stream.is_open()) out_file_stream.close();

    return 0;
}
