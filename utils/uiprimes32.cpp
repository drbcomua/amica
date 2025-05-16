#include <cstdint>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

int main() {
    constexpr uint32_t MAX_PRIME = 0xFFFFFFFF; // 2^32 - 1
    constexpr uint32_t SIEVE_LIMIT = (MAX_PRIME - 3) / 2 + 1; // Count of odd numbers from 3 to MAX_PRIME

    std::vector is_prime(SIEVE_LIMIT, true); // is_prime[i] -> 2*i + 3 is prime

    const auto sqrt_limit = static_cast<uint32_t>(std::sqrt(MAX_PRIME));

    for (uint32_t i = 0; 2 * i + 3 <= sqrt_limit; ++i) {
        if (is_prime[i]) {
            const uint32_t p = 2 * i + 3;
            // Start marking from p*p
            uint32_t start = (p * p - 3) / 2;
            for (uint32_t j = start; j < SIEVE_LIMIT; j += p) {
                is_prime[j] = false;
            }
        }
    }

    std::ofstream out("uiprimes32.dat", std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open output file.\n";
        return 1;
    }

    // Write 2 separately
    uint32_t prime = 2;
    out.write(reinterpret_cast<char*>(&prime), sizeof(uint32_t));

    for (uint32_t i = 0; i < SIEVE_LIMIT; ++i) {
        if (is_prime[i]) {
            prime = 2 * i + 3;
            out.write(reinterpret_cast<char*>(&prime), sizeof(uint32_t));
        }
    }

    out.close();
    std::cout << "Prime generation complete.\n";
    return 0;
}
