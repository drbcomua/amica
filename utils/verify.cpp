#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

// 64-bit modular multiply (to avoid overflow)
uint64_t mul64(uint64_t a, uint64_t b, uint64_t m) {
    __uint128_t z = ( __uint128_t)a * b;
    return (uint64_t)(z % m);
}

// 64-bit modular exponentiation
uint64_t powmod(uint64_t a, uint64_t e, uint64_t m) {
    uint64_t res = 1;
    while (e) {
        if (e & 1) res = mul64(res, a, m);
        a = mul64(a, a, m);
        e >>= 1;
    }
    return res;
}

bool isPrime64(uint64_t n) {
    if (n < 2) return false;
    for (uint64_t p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL}) {
        if (n == p) return true;
        if (n % p == 0) return false;
    }
    // write n-1 as d*2^s
    uint64_t d = n - 1, s = 0;
    while ((d & 1) == 0) {
        d >>= 1; s++;
    }
    // these bases suffice for determinism on 64-bit n
    for (uint64_t a : {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL}) {
        if (a % n == 0) continue;
        uint64_t x = powmod(a, d, n);
        if (x==1 || x==n-1) continue;
        bool composite = true;
        for (uint64_t r = 1; r < s; ++r) {
            x = mul64(x, x, n);
            if (x == n-1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

int main() {
    const uint64_t SLICE = 0xFFFFFFFFULL;
    std::ifstream in("FF/FF/FF/FF.dat", std::ios::binary);
    if (!in) {
        std::cerr << "Cannot open slice file\n";
        return 1;
    }

    uint32_t lo32;
    size_t idx = 0;
    while (in.read(reinterpret_cast<char*>(&lo32), sizeof(lo32))) {
        uint64_t num = (SLICE << 32) | lo32;
        if (!isPrime64(num)) {
            std::cout << "Composite at index " << idx
                      << ": 0x" << std::hex << num << "\n";
            return 1;
        }
        ++idx;
    }

    std::cout << "All " << idx << " numbers are prime.\n";
    return 0;
}
