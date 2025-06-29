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
#include <cstdlib>
#include <future>
#include <functional>
#include <map>

// Metal-cpp headers
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>


//------------------------------------------------------------------------------
// Utility and Forward Declarations
//------------------------------------------------------------------------------

// Global list of 32-bit primes
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

// Function forward declarations
std::vector<std::pair<uint64_t, uint32_t>> factor(uint64_t n_val);
std::string format_factors(const std::vector<std::pair<uint64_t, uint32_t>>& factors);
std::string classify_amicable_pair(uint64_t n_val, uint64_t s_val,
                                     const std::vector<std::pair<uint64_t, uint32_t>>& n_factors,
                                     const std::vector<std::pair<uint64_t, uint32_t>>& s_factors);
uint64_t sum_proper_divisors_cpu(uint64_t n);

//------------------------------------------------------------------------------
// Prime Generation (Optimized and Scalable)
//------------------------------------------------------------------------------
void generate_primes(uint64_t limit) {
    ::primes.clear();
    std::cout << "Generating primes up to " << limit << "..." << std::endl;
    if (limit < 2) {
        std::cout << "Prime generation complete. Found 0 primes." << std::endl;
        return;
    }
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
        if (is_prime[p] && p <= 0xFFFFFFFF) {
            ::primes.push_back(static_cast<uint32_t>(p));
        }
    }
    std::cout << "Prime generation complete. Found " << ::primes.size() << " primes." << std::endl;
}

//------------------------------------------------------------------------------
// CPU-based Helpers
//------------------------------------------------------------------------------
uint64_t sum_proper_divisors_cpu(uint64_t n) {
    if (n < 2) return 0;
    uint64_t original_n = n;
    uint64_t sum_all_divs = 1;
    uint64_t temp_n = n;
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
    if (temp_n > 1) sum_all_divs *= (1 + temp_n);
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
    if (temp_n > 1) factors_list.emplace_back(temp_n, 1);
    return factors_list;
}

std::vector<std::pair<uint64_t, uint32_t>> derive_gcd_factors(const std::vector<std::pair<uint64_t, uint32_t>>& f1, const std::vector<std::pair<uint64_t, uint32_t>>& f2) {
    std::vector<std::pair<uint64_t, uint32_t>> gcd_f;
    auto it1 = f1.begin(); auto it2 = f2.begin();
    while (it1 != f1.end() && it2 != f2.end()) {
        if (it1->first < it2->first) ++it1;
        else if (it2->first < it1->first) ++it2;
        else { gcd_f.emplace_back(it1->first, std::min(it1->second, it2->second)); ++it1; ++it2; }
    }
    return gcd_f;
}

std::vector<std::pair<uint64_t, uint32_t>> derive_quotient_factors(const std::vector<std::pair<uint64_t, uint32_t>>& num_f, const std::vector<std::pair<uint64_t, uint32_t>>& div_f) {
    std::vector<std::pair<uint64_t, uint32_t>> q_f;
    auto it_n = num_f.begin(); auto it_d = div_f.begin();
    while (it_n != num_f.end()) {
        if (it_d == div_f.end() || it_n->first < it_d->first) { q_f.push_back(*it_n); ++it_n; }
        else if (it_n->first > it_d->first) ++it_d;
        else { if (it_n->second > it_d->second) q_f.emplace_back(it_n->first, it_n->second - it_d->second); ++it_n; ++it_d; }
    }
    return q_f;
}

std::string classify_amicable_pair(uint64_t n, uint64_t s, const std::vector<std::pair<uint64_t, uint32_t>>& n_f, const std::vector<std::pair<uint64_t, uint32_t>>& s_f) {
    uint64_t g = std::gcd(n, s);
    auto g_f = derive_gcd_factors(n_f, s_f);
    auto M_f = derive_quotient_factors(n_f, g_f);
    auto N_f = derive_quotient_factors(s_f, g_f);
    bool sqfM = true; for (const auto& f : M_f) if (f.second > 1) { sqfM = false; break; }
    bool sqfN = true; for (const auto& f : N_f) if (f.second > 1) { sqfN = false; break; }
    bool cpM = (std::gcd(n / g, g) == 1); bool cpN = (std::gcd(s / g, g) == 1);
    return (sqfM && sqfN && cpM && cpN ? "" : "X") + std::to_string(M_f.size()) + "," + std::to_string(N_f.size());
}

std::string format_factors(const std::vector<std::pair<uint64_t, uint32_t>>& factors) {
    std::ostringstream oss;
    for (size_t i = 0; i < factors.size(); ++i) {
        oss << factors[i].first;
        if (factors[i].second > 1) oss << '^' << factors[i].second;
        if (i + 1 < factors.size()) oss << '*';
    }
    return oss.str();
}

//------------------------------------------------------------------------------
// Main Program Logic
//------------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <max_n (uint64)>\n";
        return 1;
    }
    uint64_t max_n_arg = 0;
    try {
        max_n_arg = std::stoull(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid max_n value '" << argv[1] << "'. " << e.what() << "\n"; return 1;
    }
    if (max_n_arg < 220) {
        std::cout << "max_n (" << max_n_arg << ") is too small. Minimum is 220.\n"; return 0;
    }

    uint64_t prime_limit = static_cast<uint64_t>(std::sqrt(static_cast<double>(max_n_arg))) + 1;
    generate_primes(prime_limit);

    NS::AutoreleasePool* pAutoreleasePool = NS::AutoreleasePool::alloc()->init();

    MTL::Device* pDevice = MTL::CreateSystemDefaultDevice();
    if (!pDevice) { std::cerr << "Error: Could not create Metal device.\n"; return 1; }
    std::cout << "Using GPU: " << pDevice->name()->utf8String() << std::endl;
    MTL::CommandQueue* pCommandQueue = pDevice->newCommandQueue();

    std::ifstream kernel_file("amicable_kernel.metal");
    if (!kernel_file.is_open()) { std::cerr << "Error: Could not open kernel file amicable_kernel.metal\n"; return 1; }
    std::string source_str(std::istreambuf_iterator<char>(kernel_file), {});

    NS::Error* pError = nullptr;
    MTL::Library* pLibrary = pDevice->newLibrary(NS::String::string(source_str.c_str(), NS::UTF8StringEncoding), nullptr, &pError);
    if (!pLibrary) { std::cerr << "Kernel library creation error: " << pError->localizedDescription()->utf8String() << std::endl; return 1; }

    MTL::Function* pFunction = pLibrary->newFunction(NS::String::string("sum_divisors_kernel", NS::UTF8StringEncoding));
    MTL::ComputePipelineState* pPipelineState = pDevice->newComputePipelineState(pFunction, &pError);
    if (!pPipelineState) { std::cerr << "Compute pipeline creation error: " << pError->localizedDescription()->utf8String() << std::endl; return 1; }

    std::cout << "Transferring " << ::primes.size() << " primes to the device.\n";
    MTL::Buffer* pPrimesBuffer = pDevice->newBuffer(::primes.data(), ::primes.size() * sizeof(uint32_t), MTL::ResourceStorageModeShared);
    uint32_t num_gpu_primes = ::primes.size();

    std::vector<AmicablePairOutput> all_found_pairs;
    const size_t batch_size = 8 * 1024 * 1024;
    std::cout << "GPU Processing in batches of " << batch_size << " elements." << std::endl;

    for (uint64_t batch_start = 1; batch_start <= max_n_arg; batch_start += batch_size) {
        uint64_t current_batch_end = std::min(batch_start + batch_size - 1, max_n_arg);
        size_t current_batch_size = static_cast<size_t>(current_batch_end - batch_start + 1);
        if (current_batch_size == 0) continue;
        std::cout << "\rProcessing batch: " << batch_start << " to " << current_batch_end << "..." << std::flush;

        MTL::CommandBuffer* pCommandBuffer = pCommandQueue->commandBuffer();
        MTL::ComputeCommandEncoder* pComputeEncoder = pCommandBuffer->computeCommandEncoder();
        pComputeEncoder->setComputePipelineState(pPipelineState);

        MTL::Buffer* pNumbersBuffer = pDevice->newBuffer(current_batch_size * sizeof(uint64_t), MTL::ResourceStorageModeShared);
        MTL::Buffer* pResultsBuffer = pDevice->newBuffer(current_batch_size * sizeof(uint64_t), MTL::ResourceStorageModeShared);
        uint64_t* n_ptr = static_cast<uint64_t*>(pNumbersBuffer->contents());
        std::iota(n_ptr, n_ptr + current_batch_size, batch_start);

        pComputeEncoder->setBuffer(pNumbersBuffer, 0, 0);
        pComputeEncoder->setBuffer(pResultsBuffer, 0, 1);
        pComputeEncoder->setBuffer(pPrimesBuffer, 0, 2);
        pComputeEncoder->setBytes(&num_gpu_primes, sizeof(uint32_t), 3);

        MTL::Size gridSize = MTL::Size(current_batch_size, 1, 1);
        NS::UInteger max_threads = pPipelineState->maxTotalThreadsPerThreadgroup();
        NS::UInteger threadGroupSize = (max_threads > current_batch_size) ? current_batch_size : max_threads;
        MTL::Size threadgroupSize = MTL::Size(threadGroupSize, 1, 1);

        pComputeEncoder->dispatchThreads(gridSize, threadgroupSize);
        pComputeEncoder->endEncoding();
        pCommandBuffer->commit();
        pCommandBuffer->waitUntilCompleted();

        uint64_t* sums_from_batch_ptr = static_cast<uint64_t*>(pResultsBuffer->contents());
        for (size_t i = 0; i < current_batch_size; ++i) {
            uint64_t n = n_ptr[i];
            uint64_t m = sums_from_batch_ptr[i];
            if (m > n) {
                uint64_t sn_of_m = sum_proper_divisors_cpu(m);
                if (sn_of_m == n) {
                    auto nf = factor(n);
                    auto sf = factor(m);
                    auto cls = classify_amicable_pair(n, m, nf, sf);
                    all_found_pairs.push_back({n, m, cls, format_factors(nf), format_factors(sf)});
                }
            }
        }

        pNumbersBuffer->release();
        pResultsBuffer->release();
    }
    std::cout << "\nAll batches processed." << std::endl;

    std::sort(all_found_pairs.begin(), all_found_pairs.end());
    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) { std::cerr << "Error: cannot open output file amicable_pairs.txt\n"; }

    std::cout << "\n--- Found " << all_found_pairs.size() << " Amicable Pairs ---\n\n";
    for (const auto& p : all_found_pairs) {
        std::ostringstream oss;
        oss << p.classification_str << '\n' << p.n << '=' << p.n_factors_str << '\n' << p.s << '=' << p.s_factors_str << "\n\n";
        std::cout << oss.str();
        if (out_file_stream.is_open()) { out_file_stream << oss.str(); }
    }
    std::cout << "Done. Found a total of " << all_found_pairs.size() << " pairs. Results also saved in amicable_pairs.txt\n";
    if (out_file_stream.is_open()) out_file_stream.close();

    pPrimesBuffer->release();
    pFunction->release();
    pLibrary->release();
    pPipelineState->release();
    pCommandQueue->release();
    pDevice->release();
    pAutoreleasePool->release();

    return 0;
}