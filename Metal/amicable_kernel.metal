#include <metal_stdlib>
using namespace metal;

/**
 * @file amicable_kernel.metal
 * @brief Metal kernel for calculating the sum of proper divisors.
 *
 * This kernel is written in Metal Shading Language (MSL) and is executed
 * on the Apple Silicon GPU. Each thread processes one number.
 */

/**
 * @brief Calculates the sum of proper divisors for a batch of numbers in parallel.
 *
 * This function is the core of the GPU computation. Each instance of this kernel
 * is executed by a separate GPU thread and processes one number from the input buffer.
 *
 * @param numbers       A read-only buffer of 64-bit unsigned integers (n). Because of
 * Unified Memory, this points to the same memory the CPU uses.
 * @param results       A write-only buffer to store the calculated sums (s). This also
 * points to shared memory.
 * @param primes        A read-only buffer containing pre-calculated prime numbers needed
 * for factorization.
 * @param num_primes    The total number of primes in the primes buffer.
 * @param gid           The globally unique ID for this thread, provided by the Metal runtime.
 * It's used as an index to determine which number to process.
 */
kernel void sum_divisors_kernel(device const ulong* numbers [[buffer(0)]],
                                device ulong* results       [[buffer(1)]],
                                device const uint* primes   [[buffer(2)]],
                                constant uint& num_primes   [[buffer(3)]],
                                uint gid [[thread_position_in_grid]]) {

    // Fetch the number 'n' for this thread to process.
    ulong n = numbers[gid];

    // The sum of proper divisors for numbers less than 2 is 0.
    if (n < 2) {
        results[gid] = 0;
        return;
    }

    ulong original_n = n;
    ulong sum_all_divs = 1; // Start with 1 as 1 is a divisor of every number.
    ulong temp_n = n;       // A temporary copy of n that will be factored down.

    // Iterate through the pre-calculated primes to find the prime factors of n.
    for (uint i = 0; i < num_primes; ++i) {
        ulong p = primes[i]; // Current prime from the global list

        // Optimization: If p*p > temp_n, the remaining temp_n must either be 1 or a prime
        // number itself. This check avoids iterating through the entire primes list.
        if (p * p > temp_n) {
            break;
        }

        if (temp_n % p == 0) {
            // This block executes if p is a factor of temp_n.
            // We calculate the sum of the geometric series (1 + p + p^2 + ... + p^k).
            ulong term_sum = 1;
            ulong p_power = 1;

            do {
                temp_n /= p;
                p_power *= p;
                term_sum += p_power;
            } while (temp_n % p == 0);

            // The sum of divisors function is multiplicative, so we multiply the results
            // from each prime factor's series.
            sum_all_divs *= term_sum;
        }
    }

    // After the loop, if temp_n is still greater than 1, it means the remaining value
    // is a prime factor itself.
    if (temp_n > 1) {
        sum_all_divs *= (1 + temp_n);
    }

    // The final result is the sum of all divisors minus the number itself.
    results[gid] = sum_all_divs - original_n;
}
