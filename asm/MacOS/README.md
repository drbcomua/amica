# amica in assembly language (macOS / Apple Silicon)

`amica_fast_mac` is a hand-written AArch64 assembly port of
[`asm/Linux/amica_fast_lnx.asm`](../Linux/amica_fast_lnx.asm), targeting
macOS on Apple Silicon (developed and tested on an Apple M1, 8 cores).
It uses the same two-wave algorithm as the Linux/x86-64 version:

1. **Wave 1** - a multi-threaded, segmented, odds-only sieve of
   Eratosthenes over the full 32-bit range, counting primes and
   tracking the largest one found.
2. **Wave 2** - a multi-threaded segmented divisor sieve (DIV-free via
   Granlund-Montgomery magic-number reciprocals) that computes S(N)
   (sum of proper divisors) for every N in L2-sized windows, then scans
   for amicable pairs with an instant in-memory S(M) lookup when the
   partner falls in the same chunk, falling back to trial division
   otherwise.

## Port notes vs. the Linux/x86-64 version

- **AVX2 -> NEON.** The AVX2 popcount (lookup-table + `pshufb`) is
  replaced by NEON `CNT` + pairwise-widening adds + `UADDLV`, which is
  both simpler and needs no lookup table at all.
- **Magic-number division.** The Newton-Raphson reciprocal construction
  ports directly (`imul` -> `mul`, truncating 64-bit multiply on both
  ISAs). The one true 128/64 division needed to build the tables
  (`(2^64-1)/p`) has a zero high dividend word, so it's just a plain
  64-bit `UDIV` on AArch64 - no software 128-bit division required.
- **Bit-indexed memory ops.** x86 `BT`/`BTR` (bit-test/-reset directly
  on a memory operand) have no AArch64 equivalent; the sieve inner
  loops compute the byte offset and bit mask explicitly and use
  `BIC`/`LDRB`/`STRB` instead.
- **Output.** Apple's arm64 ABI requires variadic arguments (as used by
  `printf`/`sprintf`) to be passed packed tightly on the stack by their
  natural size and alignment - a different, easier-to-get-subtly-wrong
  scheme than the uniform 8-byte SysV slots the Linux version relies
  on. Rather than hand-encode that per call site, this port sidesteps
  libc formatting entirely: all output is built with hand-written
  decimal formatting and emitted with a single `write(2)` call per
  printed block, under the same spinlock the Linux version uses around
  its `printf` call. The textual output format is otherwise identical.
- **Atomics.** The work-stealing task queue, prime/pair counters, and
  the largest-prime max-update use LSE atomics (`LDADDAL`, `CASAL`),
  which Apple Silicon supports natively - the direct equivalent of the
  x86 `lock xadd` / `lock cmpxchg` used in the Linux version.
- **Overflow checks omitted.** The Linux version defensively guards
  every multiply/add in `SumProperDivisors` against 64-bit overflow.
  For N < 2^32 the abundancy index is bounded well within 2^64
  (sigma(N) never gets close to a 64-bit overflow), so those paths are
  provably unreachable for this program's input domain and were
  dropped to keep the port simpler.
- **Thread count.** `num_threads` is 8, matching the Apple M1's core
  count (4 performance + 4 efficiency), vs. 16 on the reference
  16-thread Intel machine. Buffer sizes scale with `num_threads`
  accordingly.

## Build

```shell
# build for macOS, Apple Silicon (arm64):
$ clang -O2 -o amica_fast_mac amica_fast_mac.s
```

No external assembler is needed - `clang`'s integrated assembler
handles AArch64 `.s` files directly, and it links against `libSystem`
(which provides `pthread_create`/`pthread_join`, `write`, and `exit`)
automatically, so no `-lpthread` or other flags are required.

## Run

```shell
$ ./amica_fast_mac
```

This finds all amicable pairs less than 2^32-1, printing each pair in
the same format as the Linux version: an Erdos `i,j` classification
line (prefixed with `X` for an irregular pair, i.e. one where the
power of 2 differs between the two members), followed by each
member's prime factorization.

## Expected result

Tested on an Apple M1 (8 cores: 4 performance + 4 efficiency), 8GB RAM.
Every count below (prime count, largest prime, and all 1043 pairs) is an
exact match with the Linux/x86-64 reference run:

```shell
$ time ./amica_fast_mac
Primes found: 203280221, Largest: 4294967291
2,1
220=2^2*5*11
284=2^2*71

X2,2
1184=2^5*37
1210=2*5*11^2

...
4,4
4280119305=3*5*7*11*59*107*587
4498673655=3*5*7*31*53*89*293

4,3
4281566032=2^4*29*97*251*379
4446000368=2^4*113*239*10289

Checked pairs in the range 2..4294967295; found: 1043 amicable pairs
./amica_fast_mac  111.27s user 1.26s system 597% cpu 18.848 total
```

(The 8-thread M1 finishes the search in ~19s wall-clock. Most of the
speedup comes from number-theoretic early aborts in `SumProperDivisors`
- a match requires the running sigma product to divide N+M exactly, and
the remaining quotient to be at least cofactor+1 - which kill most
trial-division calls before their prime scan starts. A deterministic
Miller-Rabin test (bases 2,3,5,7,11 in Montgomery arithmetic) then
settles any cofactor that survives 128 primes' worth of scanning in one
shot, instead of scanning on to sqrt(cofactor).)
