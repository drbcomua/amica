/**
 * @file amicable_kernel.cl
 * @brief OpenCL kernel for calculating the sum of proper divisors.
 *
 * This kernel is executed on the GPU. Each thread processes one number.
 */

// A simple pragma to enable 64-bit integer support, which is standard on modern GPUs.
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/**
 * @brief Calculates the sum of proper divisors for a batch of numbers in parallel.
 *
 * This function is the core of the GPU computation. Each instance of this kernel (a "work-item")
 * is executed by a separate GPU thread and processes one number from the input buffer.
 *
 * @param numbers       A read-only buffer of 64-bit unsigned integers (n) for which to
 * calculate the sum of proper divisors.
 * @param results       A write-only buffer to store the calculated sum of proper divisors (s)
 * for each corresponding number in the input buffer.
 * @param primes        A read-only buffer containing pre-calculated prime numbers, used for
 * efficient factorization. This buffer is shared by all work-items.
 * @param num_primes    The total number of primes in the primes buffer.
 */
__kernel void sum_divisors_kernel(
    __global const ulong* numbers,
    __global ulong* results,
    __global const uint* primes,
    const uint num_primes) {

    // Get the globally unique ID for this work-item (GPU thread). This ID acts as an
    // index to determine which number this specific thread should process.
    int gid = get_global_id(0);

    // Fetch the number 'n' from the input buffer for this thread to process.
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
    // This is a direct translation of the original C++ sum_proper_divisors logic.
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

            // Divide temp_n by p until it's no longer divisible.
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
