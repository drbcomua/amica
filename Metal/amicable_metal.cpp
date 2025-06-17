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

// Metal-cpp headers. These provide a C++ interface to Apple's Metal framework.
// These headers wrap the underlying Objective-C objects.
#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>


//------------------------------------------------------------------------------
// Utility and Forward Declarations
//------------------------------------------------------------------------------

// Portability wrappers for builtins for non-x86 architectures
#if defined(__GNUC__) || defined(__clang__)
#define count_trailing_zeros_u64 __builtin_ctzll
#else
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
#endif

// Global list of 32-bit primes, generated on CPU, used by GPU
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

struct ClassificationResult {
    std::string type_str;
};

// Function forward declarations for CPU-based logic
std::vector<std::pair<uint64_t, uint32_t>> factor(uint64_t n_val);
std::string format_factors(const std::vector<std::pair<uint64_t, uint32_t>>& factors);
ClassificationResult classify_amicable_pair(
    uint64_t n_val, uint64_t s_val,
    const std::vector<std::pair<uint64_t, uint32_t>>& n_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& s_factors);

//------------------------------------------------------------------------------
// Prime Generation Sieve (Multi-Threaded CPU)
//------------------------------------------------------------------------------
// This version is a portable C++ implementation without x86 SIMD specifics.
void collect_primes(const uint64_t* sieve_data, uint64_t words,
                         uint64_t low_val_in_segment,
                         std::vector<uint32_t>& outvec) {
    for (uint64_t w = 0; w < words; ++w) {
        uint64_t word = sieve_data[w];
        while (word) {
            int bit = count_trailing_zeros_u64(word);
            uint64_t idx_bit = w * 64 + bit;
            uint64_t prime_candidate = low_val_in_segment + 2 * idx_bit;
            if (prime_candidate <= 0xFFFFFFFFULL) {
                outvec.push_back(static_cast<uint32_t>(prime_candidate));
            }
            word &= word - 1; // Clear the least significant bit
        }
    }
}

// The multi-threaded prime generation logic remains the same, but calls the portable `collect_primes`.
std::vector<uint32_t> process_segment_task(
    uint32_t segment_idx_for_task,
    uint64_t total_sieve_slots_count_const,
    uint64_t segment_bits_capacity_const,
    const std::vector<uint32_t>& base_primes_const_ref) {

    const uint64_t current_segment_low_bit_global_idx = static_cast<uint64_t>(segment_idx_for_task) * segment_bits_capacity_const;
    if (current_segment_low_bit_global_idx >= total_sieve_slots_count_const) {
        return {};
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
        if ((start_multiple & 1) == 0) start_multiple += p;
        if (start_multiple > high_val_in_segment) continue;
        uint64_t bit_index = (start_multiple - low_val_in_segment) / 2;
        for (uint64_t b = bit_index; b < current_segment_active_bits; b += p) {
            sieve_local[b >> 6] &= ~(1ULL << (b & 63));
        }
    }

    double log_low_val = std::log(static_cast<double>(low_val_in_segment > 1 ? low_val_in_segment : 2.0));
    if (log_low_val < 0.1) log_low_val = 0.1;
    segment_primes_local.reserve(static_cast<size_t>(static_cast<double>(current_segment_active_bits) / log_low_val * 1.2) + 100);
    collect_primes(sieve_local.data(), words_in_segment, low_val_in_segment, segment_primes_local);
    return segment_primes_local;
}

void generate_primes_multi_threaded(int mb_segment_size_param) {
    ::primes.clear();
    std::ios_base::sync_with_stdio(false); std::cin.tie(NULL);
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
    const uint64_t total_sieve_slots_count = (MAX_SIEVE_N - 1) / 2;
    const auto num_segments = static_cast<uint32_t>((total_sieve_slots_count + segment_bits_capacity - 1) / segment_bits_capacity);
    ::primes.push_back(2);
    double log_max_n = std::log(static_cast<double>(MAX_SIEVE_N > 1 ? MAX_SIEVE_N : 2.0));
    if (log_max_n < 0.1) log_max_n = 0.1;
    ::primes.reserve(static_cast<size_t>(static_cast<double>(MAX_SIEVE_N) / log_max_n * 1.1) + 100);
    unsigned int num_hw_threads = std::thread::hardware_concurrency();
    if (num_hw_threads == 0) num_hw_threads = 1;
    num_hw_threads = std::min(num_hw_threads, static_cast<unsigned int>(32));
    num_hw_threads = std::min(num_hw_threads, num_segments > 0 ? num_segments : 1U);
    std::cout << "Generating primes in memory up to " << MAX_SIEVE_N << " (using " << mb_segment_size_param << "MB segments, " << num_hw_threads << " threads)...\n";
    std::vector<std::vector<uint32_t>> results_from_all_segments(num_segments);
    for (uint32_t batch_start_segment_idx = 0; batch_start_segment_idx < num_segments; batch_start_segment_idx += num_hw_threads) {
        std::vector<std::future<std::vector<uint32_t>>> futures_in_batch;
        futures_in_batch.reserve(num_hw_threads);
        for (unsigned int i = 0; i < num_hw_threads; ++i) {
            uint32_t current_segment_to_process = batch_start_segment_idx + i;
            if (current_segment_to_process >= num_segments) break;
            futures_in_batch.push_back(std::async(std::launch::async, process_segment_task, current_segment_to_process, total_sieve_slots_count, segment_bits_capacity, std::cref(base_primes)));
        }
        for (size_t i = 0; i < futures_in_batch.size(); ++i) {
            uint32_t actual_segment_idx_processed = batch_start_segment_idx + static_cast<uint32_t>(i);
            results_from_all_segments[actual_segment_idx_processed] = futures_in_batch[i].get();
            if (((actual_segment_idx_processed + 1) % (num_segments / 100 + 1) == 0) || (actual_segment_idx_processed == num_segments - 1)) {
                std::cout << "\rProcessed segment " << actual_segment_idx_processed + 1 << "/" << num_segments << std::flush;
            }
        }
    }
    std::cout << "\nAll segments processed. Consolidating primes...\n";
    for (const auto& segment_primes_vec : results_from_all_segments) {
        ::primes.insert(::primes.end(), segment_primes_vec.begin(), segment_primes_vec.end());
    }
    std::cout << "Prime generation complete. Found " << ::primes.size() << " primes.\n";
}


//------------------------------------------------------------------------------
// Amicable Pair Logic (CPU-based helpers)
//------------------------------------------------------------------------------
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

ClassificationResult classify_amicable_pair(uint64_t n, uint64_t s, const std::vector<std::pair<uint64_t, uint32_t>>& n_f, const std::vector<std::pair<uint64_t, uint32_t>>& s_f) {
    ClassificationResult res;
    uint64_t g = std::gcd(n, s);
    auto g_f = derive_gcd_factors(n_f, s_f);
    auto M_f = derive_quotient_factors(n_f, g_f);
    auto N_f = derive_quotient_factors(s_f, g_f);
    bool sqfM = true; for (const auto& f : M_f) if (f.second > 1) { sqfM = false; break; }
    bool sqfN = true; for (const auto& f : N_f) if (f.second > 1) { sqfN = false; break; }
    bool cpM = (std::gcd(n / g, g) == 1); bool cpN = (std::gcd(s / g, g) == 1);
    res.type_str = (sqfM && sqfN && cpM && cpN ? "" : "X") + std::to_string(M_f.size()) + "," + std::to_string(N_f.size());
    return res;
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
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <max_n (uint64)> [segment_size_mb (optional, int, default 16)]\n";
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
    int segment_size_mb_val = 16;
    if (argc == 3) {
        try {
            segment_size_mb_val = std::stoi(argv[2]);
            if (segment_size_mb_val < 1 || segment_size_mb_val > 2048) segment_size_mb_val = 16;
        } catch (...) { segment_size_mb_val = 16; }
    }

    // --- Phase 1: Generate Primes on CPU ---
    generate_primes_multi_threaded(segment_size_mb_val);
    if (::primes.empty()) { std::cerr << "Error: Prime generation failed.\n"; return 1; }
    const uint64_t search_limit_val = max_n_arg;

    // --- Phase 2: Amicable Pair Search on GPU using Metal ---
    NS::AutoreleasePool* pAutoreleasePool = NS::AutoreleasePool::alloc()->init();

    // 1. Metal Setup
    MTL::Device* pDevice = MTL::CreateSystemDefaultDevice();
    std::cout << "Using GPU: " << pDevice->name()->utf8String() << std::endl;
    MTL::CommandQueue* pCommandQueue = pDevice->newCommandQueue();

    // 2. Load and Compile Kernel from source file
    std::ifstream kernel_file("amicable_kernel.metal");
    if (!kernel_file.is_open()) { std::cerr << "Error: Could not open kernel file amicable_kernel.metal\n"; return 1; }
    std::string source_str(std::istreambuf_iterator<char>(kernel_file), {});

    NS::Error* pError = nullptr;
    MTL::Library* pLibrary = pDevice->newLibrary(NS::String::string(source_str.c_str(), NS::UTF8StringEncoding), nullptr, &pError);
    if (!pLibrary) { std::cerr << "Kernel library creation error: " << pError->localizedDescription()->utf8String() << std::endl; return 1; }

    MTL::Function* pFunction = pLibrary->newFunction(NS::String::string("sum_divisors_kernel", NS::UTF8StringEncoding));
    MTL::ComputePipelineState* pPipelineState = pDevice->newComputePipelineState(pFunction, &pError);
    if (!pPipelineState) { std::cerr << "Compute pipeline creation error: " << pError->localizedDescription()->utf8String() << std::endl; return 1; }

    // 3. Create Shared Buffers
    uint64_t sqrt_search_limit = static_cast<uint64_t>(std::sqrt(static_cast<double>(search_limit_val))) + 1;
    std::vector<uint32_t> gpu_primes;
    gpu_primes.reserve(static_cast<size_t>(sqrt_search_limit / std::log(sqrt_search_limit > 1 ? sqrt_search_limit : 2.0)));
    for (uint32_t p : ::primes) {
        if (p > sqrt_search_limit) break;
        gpu_primes.push_back(p);
    }
    std::cout << "Optimizing for Metal: Using " << gpu_primes.size() << " primes (up to sqrt(" << search_limit_val << ")) in shared memory.\n";
    uint32_t num_gpu_primes = gpu_primes.size();

    // Create the results buffer in SHARED memory, accessible by both CPU and GPU
    MTL::Buffer* pResultsBuffer = pDevice->newBuffer((search_limit_val + 1) * sizeof(uint64_t), MTL::ResourceStorageModeShared);
    uint64_t* s_results_ptr = static_cast<uint64_t*>(pResultsBuffer->contents());

    // --- GPU Computation in Batches ---
    const size_t batch_size = 2 * 1024 * 1024; // 2 million elements per batch is a safe size
    std::cout << "GPU Processing in batches of " << batch_size << " elements." << std::endl;

    for (uint64_t batch_start = 0; batch_start <= search_limit_val; batch_start += batch_size) {
        uint64_t current_batch_end = std::min(batch_start + batch_size - 1, search_limit_val);
        size_t current_batch_size = static_cast<size_t>(current_batch_end - batch_start + 1);
        if (current_batch_size == 0) continue;
        std::cout << "\rProcessing batch: " << batch_start << " to " << current_batch_end << "..." << std::flush;

        MTL::CommandBuffer* pCommandBuffer = pCommandQueue->commandBuffer();
        MTL::ComputeCommandEncoder* pComputeEncoder = pCommandBuffer->computeCommandEncoder();

        pComputeEncoder->setComputePipelineState(pPipelineState);

        // Create buffers for this batch
        MTL::Buffer* pNumbersBuffer = pDevice->newBuffer(current_batch_size * sizeof(uint64_t), MTL::ResourceStorageModeShared);
        uint64_t* n_ptr = static_cast<uint64_t*>(pNumbersBuffer->contents());
        std::iota(n_ptr, n_ptr + current_batch_size, batch_start);

        MTL::Buffer* pPrimesBuffer = pDevice->newBuffer(gpu_primes.data(), gpu_primes.size() * sizeof(uint32_t), MTL::ResourceStorageModeShared);

        pComputeEncoder->setBuffer(pNumbersBuffer, 0, 0);
        pComputeEncoder->setBuffer(pResultsBuffer, (batch_start * sizeof(uint64_t)), 1); // IMPORTANT: Set correct offset for results
        pComputeEncoder->setBuffer(pPrimesBuffer, 0, 2);
        pComputeEncoder->setBytes(&num_gpu_primes, sizeof(uint32_t), 3);

        MTL::Size gridSize = MTL::Size(current_batch_size, 1, 1);
        NS::UInteger threadGroupSize = pPipelineState->maxTotalThreadsPerThreadgroup();
        if (threadGroupSize > current_batch_size) threadGroupSize = current_batch_size;
        MTL::Size threadgroupSize = MTL::Size(threadGroupSize, 1, 1);

        pComputeEncoder->dispatchThreads(gridSize, threadgroupSize);
        pComputeEncoder->endEncoding();

        pCommandBuffer->commit();
        pCommandBuffer->waitUntilCompleted(); // Wait for this batch to finish

        pNumbersBuffer->release();
        pPrimesBuffer->release();
    }
    std::cout << "\nAll batches processed." << std::endl;

    // --- Filter candidates on CPU using the shared results buffer ---
    std::cout << "Filtering potential pairs on CPU..." << std::endl;
    std::vector<AmicablePairOutput> all_found_pairs;
    for (uint64_t n = 220; n <= search_limit_val; ++n) {
        uint64_t s = s_results_ptr[n];
        if (s > n && s <= search_limit_val) {
            if (s_results_ptr[s] == n) {
                auto nf = factor(n);
                auto sf = factor(s);
                auto cls = classify_amicable_pair(n, s, nf, sf);
                all_found_pairs.push_back({n, s, cls.type_str, format_factors(nf), format_factors(sf)});
            }
        }
    }
    std::cout << "GPU computation and filtering complete." << std::endl;

    // --- Sort and Print Results ---
    std::sort(all_found_pairs.begin(), all_found_pairs.end());
    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) { std::cerr << "Error: cannot open output file amicable_pairs.txt\n"; }
    for (const auto& p : all_found_pairs) {
        std::ostringstream oss;
        oss << p.classification_str << '\n' << p.n << '=' << p.n_factors_str << '\n' << p.s << '=' << p.s_factors_str << "\n\n";
        std::cout << oss.str();
        if (out_file_stream.is_open()) { out_file_stream << oss.str(); }
    }
    std::cout << "Done. Found " << all_found_pairs.size() << " pairs. Results in amicable_pairs.txt\n";
    if (out_file_stream.is_open()) out_file_stream.close();

    // --- Cleanup Metal objects ---
    pResultsBuffer->release();
    pFunction->release();
    pLibrary->release();
    pPipelineState->release();
    pCommandQueue->release();
    pDevice->release();
    pAutoreleasePool->release();

    return 0;
}
