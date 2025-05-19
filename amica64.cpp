#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <numeric>

//------------------------------------------------------------------------------
// Global list of 32-bit primes
static std::vector<uint32_t> primes;

// Load primes from your binary file
void load_primes(const std::string& filename) {
    std::ifstream in{filename, std::ios::binary};
    if (!in) {
        std::cerr << "Error: cannot open primes file '"
                  << filename << "'\n";
        std::exit(1);
    }
    uint32_t p;
    while (in.read(reinterpret_cast<char*>(&p), sizeof(p))) {
        primes.push_back(p);
    }
    std::cout << "Loaded " << primes.size()
              << " primes from " << filename << "\n";
}

// Sum of proper divisors via σ(n) formula
uint64_t sum_proper_divisors(uint64_t n) {
    if (n < 2) return 0;
    const uint64_t original = n;
    uint64_t sigma = 1;
    for (uint32_t p : primes) {
        uint64_t pp = p;
        if (pp*pp > n) break;
        if (n % pp == 0) {
            uint64_t term = 1, power = 1;
            while (n % pp == 0) {
                n /= pp;
                power *= pp;
                term += power;
            }
            sigma *= term;
        }
    }
    if (n > 1) sigma *= 1 + n;
    return sigma - original;
}

// Factor x into (prime, exponent) list
std::vector<std::pair<uint64_t,uint32_t>> factor(uint64_t x) {
    std::vector<std::pair<uint64_t,uint32_t>> factors;
    uint64_t n = x;
    for (uint32_t p : primes) {
        uint64_t pp = p;
        if (pp*pp > n) break;
        if (n % pp == 0) {
            uint32_t exp = 0;
            while (n % pp == 0) {
                n /= pp;
                ++exp;
            }
            factors.emplace_back(pp, exp);
        }
    }
    if (n > 1) factors.emplace_back(n,1);
    return factors;
}

// Format factor list as "p1^e1*p2^e2*…"
std::string format_factors(const std::vector<std::pair<uint64_t,uint32_t>>& factors) {
    std::string s;
    for (size_t i = 0; i < factors.size(); ++i) {
        auto [p,e] = factors[i];
        s += std::to_string(p);
        if (e > 1) {
            s += '^';
            s += std::to_string(e);
        }
        if (i+1 < factors.size()) s += '*';
    }
    return s;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <max_n (uint64)> [primes_file]\n";
        return 1;
    }
    uint64_t max_n = std::stoull(argv[1]);
    std::string primes_file = argc>2 ? argv[2] : "uiprimes32.dat";

    load_primes(primes_file);

    std::ofstream out{"amicable_pairs.txt"};
    if (!out) {
        std::cerr << "Error: cannot open output file\n";
        return 1;
    }

    for (uint64_t n = 2; n <= max_n; ++n) {
        uint64_t s = sum_proper_divisors(n);
        if (s > n && s <= max_n && sum_proper_divisors(s) == n) {
            // 1) classification
            uint64_t g = std::gcd(n, s);
            uint64_t M = n / g, N = s / g;
            auto fM = factor(M), fN = factor(N);

            bool sqfM = true, sqfN = true;
            for (auto &pr : fM) if (pr.second > 1) { sqfM = false; break; }
            for (auto &pr : fN) if (pr.second > 1) { sqfN = false; break; }

            bool cpM = (std::gcd(M, g) == 1);
            bool cpN = (std::gcd(N, g) == 1);

            bool regular = (sqfM && sqfN && cpM && cpN);
            std::string prefix = regular ? "" : "X";
            size_t i = fM.size(), j = fN.size();

            // 2) print classification line
            out   << prefix << i << ',' << j << '\n';
            std::cout << prefix << i << ',' << j << '\n';

            // 3) print the two factorizations
            auto fn = factor(n), fs = factor(s);
            out   << n << '=' << format_factors(fn) << '\n';
            out   << s << '=' << format_factors(fs) << "\n\n";

            std::cout << n << '=' << format_factors(fn) << '\n'
                      << s << '=' << format_factors(fs) << "\n\n";
        }
    }

    std::cout << "Done. Results in amicable_pairs.txt\n";
    return 0;
}
