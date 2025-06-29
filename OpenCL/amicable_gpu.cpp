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

// Define the OpenCL target version BEFORE including the header.
#define CL_TARGET_OPENCL_VERSION 200
#include <CL/cl.h>

// A simple macro to check for OpenCL errors.
#define CL_CHECK(err) \
    do { \
        if (err != CL_SUCCESS) { \
            std::cerr << "OpenCL Error: " << err << " at line " << __LINE__ << std::endl; \
            exit(1); \
        } \
    } while (0)

// Global list of primes, used by CPU factorizer and GPU kernel
static std::vector<uint32_t> primes;

// Struct to hold all data for a found pair for sorted output
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
// Prime Generation (Optimized)
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

    if (temp_n > 1) {
        sum_all_divs *= (1 + temp_n);
    }

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
        else { gcd_factors.emplace_back(it1->first, std::min(it1->second, it2->second)); ++it1; ++it2; }
    }
    return gcd_factors;
}

std::vector<std::pair<uint64_t, uint32_t>> derive_quotient_factors(
    const std::vector<std::pair<uint64_t, uint32_t>>& num_factors,
    const std::vector<std::pair<uint64_t, uint32_t>>& divisor_factors) {
    std::vector<std::pair<uint64_t, uint32_t>> quotient_factors;
    auto it_num = num_factors.begin(); auto it_div = divisor_factors.begin();
    while (it_num != num_factors.end()) {
        if (it_div == divisor_factors.end() || it_num->first < it_div->first) { quotient_factors.push_back(*it_num); ++it_num; }
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

std::string format_factors(const std::vector<std::pair<uint64_t, uint32_t>>& factors) {
    std::ostringstream oss;
    for (size_t i = 0; i < factors.size(); ++i) {
        auto [p, e] = factors[i]; oss << p;
        if (e > 1) { oss << '^' << e; }
        if (i + 1 < factors.size()) { oss << '*'; }
    }
    return oss.str();
}

//------------------------------------------------------------------------------
// OpenCL Kernel Source
//------------------------------------------------------------------------------
const char* kernel_source = R"CLC(
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
    __kernel void sum_divisors_kernel(
        __global const ulong* numbers, __global ulong* results,
        __global const uint* primes, const uint num_primes) {
        int gid = get_global_id(0);
        ulong n = numbers[gid];
        if (n < 2) { results[gid] = 0; return; }
        ulong original_n = n;
        ulong sum_all_divs = 1;
        ulong temp_n = n;
        for (uint i = 0; i < num_primes; ++i) {
            ulong p = primes[i];
            if (p * p > temp_n) break;
            if (temp_n % p == 0) {
                ulong term_sum = 1;
                ulong p_power = 1;
                do {
                    temp_n /= p;
                    p_power *= p;
                    term_sum += p_power;
                } while (temp_n % p == 0);
                sum_all_divs *= term_sum;
            }
        }
        if (temp_n > 1) { sum_all_divs *= (1 + temp_n); }
        results[gid] = sum_all_divs - original_n;
    }
)CLC";

//------------------------------------------------------------------------------
// Main Program Logic
//------------------------------------------------------------------------------
int main(const int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);

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
        std::cout << "max_n (" << max_n_arg << ") is too small. Minimum for amicable pairs is 220.\n";
        return 0;
    }

    uint64_t prime_limit = static_cast<uint64_t>(std::sqrt(static_cast<double>(max_n_arg))) + 1;
    generate_primes(prime_limit);

    cl_int err;
    cl_platform_id platform_id;
    cl_device_id device_id;

    err = clGetPlatformIDs(1, &platform_id, NULL); CL_CHECK(err);
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);
    if (err == CL_DEVICE_NOT_FOUND) {
        std::cout << "GPU not found. Trying OpenCL on CPU..." << std::endl;
        err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_CPU, 1, &device_id, NULL); CL_CHECK(err);
    }

    char deviceName[1024];
    clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(deviceName), deviceName, NULL);
    std::cout << "Using OpenCL device: " << deviceName << std::endl;

    cl_context context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &err); CL_CHECK(err);
    cl_command_queue command_queue = clCreateCommandQueueWithProperties(context, device_id, NULL, &err); CL_CHECK(err);

    cl_program program = clCreateProgramWithSource(context, 1, &kernel_source, NULL, &err); CL_CHECK(err);
    std::cout << "Building OpenCL kernel..." << std::endl;
    err = clBuildProgram(program, 1, &device_id, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t log_size;
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
        std::vector<char> log(log_size);
        clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, log_size, log.data(), NULL);
        std::cerr << "Kernel Compilation Error:\n" << log.data() << std::endl;
        return 1;
    }
    cl_kernel kernel = clCreateKernel(program, "sum_divisors_kernel", &err); CL_CHECK(err);

    std::cout << "Transferring " << ::primes.size() << " primes to the device." << std::endl;
    cl_mem primes_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, ::primes.size() * sizeof(cl_uint), ::primes.data(), &err); CL_CHECK(err);
    cl_uint num_primes_val = static_cast<cl_uint>(::primes.size());

    std::vector<AmicablePairOutput> all_found_pairs;
    const size_t BATCH_SIZE_ELEMENTS = 8 * 1024 * 1024;

    std::cout << "Processing numbers up to " << max_n_arg << " on the device..." << std::endl;
    for (uint64_t batch_start = 1; batch_start <= max_n_arg; batch_start += BATCH_SIZE_ELEMENTS) {
        uint64_t current_batch_end = std::min(batch_start + BATCH_SIZE_ELEMENTS - 1, max_n_arg);
        size_t current_batch_size = static_cast<size_t>(current_batch_end - batch_start + 1);

        std::cout << "\rProcessing batch: " << batch_start << " to " << current_batch_end << "..." << std::flush;

        std::vector<uint64_t> numbers_in_batch(current_batch_size);
        std::iota(numbers_in_batch.begin(), numbers_in_batch.end(), batch_start);
        std::vector<uint64_t> sums_from_batch(current_batch_size);

        cl_mem n_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, numbers_in_batch.size() * sizeof(cl_ulong), numbers_in_batch.data(), &err); CL_CHECK(err);
        cl_mem s_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sums_from_batch.size() * sizeof(cl_ulong), NULL, &err); CL_CHECK(err);

        clSetKernelArg(kernel, 0, sizeof(cl_mem), &n_buf);
        clSetKernelArg(kernel, 1, sizeof(cl_mem), &s_buf);
        clSetKernelArg(kernel, 2, sizeof(cl_mem), &primes_buf);
        clSetKernelArg(kernel, 3, sizeof(cl_uint), &num_primes_val);

        err = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL, &current_batch_size, NULL, 0, NULL, NULL); CL_CHECK(err);
        err = clEnqueueReadBuffer(command_queue, s_buf, CL_TRUE, 0, sums_from_batch.size() * sizeof(cl_ulong), sums_from_batch.data(), 0, NULL, NULL); CL_CHECK(err);

        clReleaseMemObject(n_buf);
        clReleaseMemObject(s_buf);

        for (size_t i = 0; i < current_batch_size; ++i) {
            uint64_t n = numbers_in_batch[i];
            uint64_t m = sums_from_batch[i];

            if (m > n) {
                uint64_t sn_of_m = sum_proper_divisors_cpu(m);
                if (sn_of_m == n) {
                    auto nf = factor(n);
                    auto sf = factor(m);
                    auto cls = classify_amicable_pair(n, m, nf, sf);
                    all_found_pairs.push_back({ n, m, cls, format_factors(nf), format_factors(sf) });
                }
            }
        }
    }
    std::cout << "\nAll batches processed." << std::endl;

    std::sort(all_found_pairs.begin(), all_found_pairs.end());

    std::ofstream out_file_stream{"amicable_pairs.txt"};
    if (!out_file_stream) {
        std::cerr << "Error: cannot open output file amicable_pairs.txt\n";
    }

    std::cout << "\n--- Found " << all_found_pairs.size() << " Amicable Pairs ---\n\n";
    for (const auto& p : all_found_pairs) {
        std::ostringstream oss;
        oss << p.classification_str << '\n' << p.n << '=' << p.n_factors_str << '\n' << p.s << '=' << p.s_factors_str << "\n\n";
        std::cout << oss.str();
        if(out_file_stream.is_open()){
            out_file_stream << oss.str();
        }
    }

    // *** MODIFIED LINE FOR FINAL SUMMARY ***
    std::cout << "Done. Found a total of " << all_found_pairs.size() << " pairs. Results also saved in amicable_pairs.txt\n";

    if (out_file_stream.is_open()) out_file_stream.close();

    clReleaseMemObject(primes_buf);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);

    return 0;
}