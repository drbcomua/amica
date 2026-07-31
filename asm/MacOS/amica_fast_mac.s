// -------------------------------------------------------------------
// amica_fast_mac.s - Multi-threaded Segmented Sieve for macOS arm64
// Supports full 32-bit unsigned range (0 to 4,294,967,295)
// Optimized for Apple M1 (8 cores: 4 performance + 4 efficiency)
// Ported from asm/Linux/amica_fast_lnx.asm (x86-64/AVX2) to AArch64.
//
// Port notes vs the Linux/x86-64 version:
//  - AVX2 popcount (lookup-table + pshufb) replaced by NEON CNT+UADDLV.
//  - Magic-number reciprocal division (Newton-Raphson) ported directly:
//    x86 MUL (low 64) -> AArch64 MUL; the one true 128/64 DIV needed
//    (computing (2^64-1)/p) has a zero high dividend word, so it is a
//    plain 64-bit UDIV on AArch64.
//  - x86 BT/BTR bit-indexed memory ops have no AArch64 equivalent;
//    replaced with explicit byte/bit address computation + BIC/ORR.
//  - libc printf/sprintf are variadic; Apple's arm64 ABI packs variadic
//    args tightly on the stack by natural size/alignment (unlike the
//    uniform 8-byte slots of SysV x86-64), which is easy to get subtly
//    wrong by hand. To keep the port robust, all output is generated
//    with hand-written decimal formatting and a single write(2) syscall
//    per printed line/block, under the same spinlock the Linux version
//    uses around its printf call.
//  - Thread/global counters use LSE atomics (LDADDAL/CASAL), which
//    Apple Silicon supports natively (equivalent to `lock xadd`/`lock
//    cmpxchg`).
// -------------------------------------------------------------------

.equ limit,          4294967295
.equ sqrt_limit,      262144      // Cover factors up to 68.7 billion for S(N)
.equ segment_size,    131072      // 128 KB (Fits in L2 cache)
.equ num_threads,     8           // Matches Apple M1 (4P+4E cores)
.equ sub_seg_elems,   32768       // Divisor-sieve window: 128KB rem + 256KB sigma
.equ chunk_size,      2097152     // 2M integers per task (matches segment_size*16 bits)
.equ max_base_primes, 30000
.equ per_thread_idx_stride, max_base_primes*8

// The AArch64 `mov reg,#imm` alias only covers immediates that fit a
// single MOVZ/MOVN/logical-immediate encoding; arbitrary 64-bit
// constants (e.g. `limit`) need an explicit MOVZ+MOVK sequence.
.macro LDIMM64 reg, val
    movz \reg, #(\val) & 0xFFFF
    movk \reg, #((\val) >> 16) & 0xFFFF, lsl #16
    movk \reg, #((\val) >> 32) & 0xFFFF, lsl #32
    movk \reg, #((\val) >> 48) & 0xFFFF, lsl #48
.endm

// =====================================================================
// __DATA,__bss  (uninitialized, zero at load)
// =====================================================================
.zerofill __DATA,__bss,_base_primes,max_base_primes*4,2
.zerofill __DATA,__bss,_base_prime_cnt,4,2
.zerofill __DATA,__bss,_sieve_prime_cnt,4,2

.zerofill __DATA,__bss,_pinv_tab,max_base_primes*8,3
.zerofill __DATA,__bss,_plim_tab,max_base_primes*8,3
.zerofill __DATA,__bss,_psq_tab,max_base_primes*8,3

.zerofill __DATA,__bss,_thread_handles,num_threads*8,3

.zerofill __DATA,__bss,_segment_buffers,num_threads*segment_size,5
.zerofill __DATA,__bss,_next_indices,num_threads*max_base_primes*8,3
.zerofill __DATA,__bss,_base_sieve_buffer,32768,5

.zerofill __DATA,__bss,_global_task_start,8,3
.zerofill __DATA,__bss,_global_prime_cnt,8,3
.zerofill __DATA,__bss,_global_pair_cnt,8,3
.zerofill __DATA,__bss,_global_largest,4,2
.zerofill __DATA,__bss,_print_lock,4,2

.zerofill __DATA,__bss,_chunk_rem,num_threads*chunk_size*4,3
.zerofill __DATA,__bss,_chunk_sigma,num_threads*chunk_size*8,3

// =====================================================================
// __TEXT,__const  (string literals)
// =====================================================================
.section __TEXT,__const
msg_primes_1: .ascii "Primes found: "
.equ msg_primes_1_len, . - msg_primes_1
msg_primes_2: .ascii ", Largest: "
.equ msg_primes_2_len, . - msg_primes_2
nl: .ascii "\n"
comma: .ascii ","
star: .ascii "*"
caret: .ascii "^"
eqsign: .ascii "="
xflag: .ascii "X"
msg_summary_1: .ascii "Checked pairs in the range 2.."
.equ msg_summary_1_len, . - msg_summary_1
msg_summary_2: .ascii "; found: "
.equ msg_summary_2_len, . - msg_summary_2
msg_summary_3: .ascii " amicable pairs\n"
.equ msg_summary_3_len, . - msg_summary_3

.align 4
iota4: .word 0, 1, 2, 3

.text
.align 2

// =====================================================================
// main
// =====================================================================
.globl _main
_main:
    stp x29, x30, [sp, #-16]!
    mov x29, sp

    bl _build_base_sieve
    bl _collect_base_primes
    bl _build_recip_tables

    // Initialize globals (bss already zero, but be explicit like the x86 version)
    adrp x0, _global_task_start@PAGE
    add  x0, x0, _global_task_start@PAGEOFF
    str  xzr, [x0]
    adrp x0, _global_prime_cnt@PAGE
    add  x0, x0, _global_prime_cnt@PAGEOFF
    str  xzr, [x0]
    adrp x0, _global_largest@PAGE
    add  x0, x0, _global_largest@PAGEOFF
    str  wzr, [x0]

    // ---------------------------------------------------------------
    // WAVE 1: Prime Sieve
    // ---------------------------------------------------------------
    bl _run_wave1

    adrp x9, _global_prime_cnt@PAGE
    add  x9, x9, _global_prime_cnt@PAGEOFF
    ldr  x0, [x9]                     // prime count
    adrp x9, _global_largest@PAGE
    add  x9, x9, _global_largest@PAGEOFF
    ldr  w1, [x9]                     // largest prime

    bl _print_wave1_summary

    // ---------------------------------------------------------------
    // WAVE 2: Amicable Pair Search (Segmented Divisor Sieve)
    // ---------------------------------------------------------------
    adrp x0, _global_task_start@PAGE
    add  x0, x0, _global_task_start@PAGEOFF
    str  xzr, [x0]

    bl _run_wave2

    adrp x9, _global_pair_cnt@PAGE
    add  x9, x9, _global_pair_cnt@PAGEOFF
    ldr  x0, [x9]
    bl _print_wave2_summary

    mov x0, #0
    bl _exit

// =====================================================================
// _build_base_sieve - sieve of Eratosthenes over odd numbers 1..sqrt_limit
// Bit k (byte k>>3, bit k&7) represents odd number n = 2k+1.
// =====================================================================
_build_base_sieve:
    stp x29, x30, [sp, #-16]!
    mov x29, sp

    adrp x9, _base_sieve_buffer@PAGE
    add  x9, x9, _base_sieve_buffer@PAGEOFF
    mov  x10, #4096
    mov  x11, #-1
Lbs_fill:
    str  x11, [x9], #8
    subs x10, x10, #1
    b.ne Lbs_fill

    adrp x9, _base_sieve_buffer@PAGE
    add  x9, x9, _base_sieve_buffer@PAGEOFF

    // n=1 (k=0) is not prime
    ldrb w12, [x9]
    and  w12, w12, #0xFE
    strb w12, [x9]

    mov  x13, #3            // candidate n
Lbs_loop:
    cmp  x13, #512
    b.hs Lbs_done

    sub  x14, x13, #1
    lsr  x14, x14, #1       // k
    lsr  x15, x14, #3
    and  x16, x14, #7
    ldrb w17, [x9, x15]
    lsr  w0, w17, w16
    and  w0, w0, #1
    cbz  w0, Lbs_next

    mul  x8, x13, x13        // start = n*n
    lsl  x1, x13, #1         // step = 2n
Lbs_clear:
    cmp  x8, #sqrt_limit
    b.hs Lbs_clear_done
    sub  x2, x8, #1
    lsr  x2, x2, #1
    lsr  x3, x2, #3
    and  x4, x2, #7
    ldrb w5, [x9, x3]
    mov  w6, #1
    lsl  w6, w6, w4
    bic  w5, w5, w6
    strb w5, [x9, x3]
    add  x8, x8, x1
    b    Lbs_clear
Lbs_clear_done:
Lbs_next:
    add  x13, x13, #2
    b    Lbs_loop
Lbs_done:
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _collect_base_primes - scans base_sieve_buffer for n=3..sqrt_limit-1,
// appends primes to base_primes[], sets base_prime_cnt and
// sieve_prime_cnt (count of primes < 65536).
// =====================================================================
_collect_base_primes:
    stp x29, x30, [sp, #-16]!
    mov x29, sp

    adrp x9, _base_sieve_buffer@PAGE
    add  x9, x9, _base_sieve_buffer@PAGEOFF
    adrp x20, _base_primes@PAGE
    add  x20, x20, _base_primes@PAGEOFF

    mov  x21, #0             // count
    mov  x13, #3              // candidate n
Lcb_loop:
    cmp  x13, #sqrt_limit
    b.hs Lcb_done

    sub  x14, x13, #1
    lsr  x14, x14, #1
    lsr  x15, x14, #3
    and  x16, x14, #7
    ldrb w17, [x9, x15]
    lsr  w0, w17, w16
    and  w0, w0, #1
    cbz  w0, Lcb_next

    str  w13, [x20, x21, lsl #2]
    add  x21, x21, #1
Lcb_next:
    add  x13, x13, #2
    b    Lcb_loop
Lcb_done:
    adrp x0, _base_prime_cnt@PAGE
    add  x0, x0, _base_prime_cnt@PAGEOFF
    str  w21, [x0]

    // sieve_prime_cnt = count of primes < 65536
    mov  x22, #0
Lcb_count:
    cmp  x22, x21
    b.hs Lcb_count_done
    ldr  w0, [x20, x22, lsl #2]
    cmp  w0, #65536
    b.hs Lcb_count_done
    add  x22, x22, #1
    b    Lcb_count
Lcb_count_done:
    adrp x0, _sieve_prime_cnt@PAGE
    add  x0, x0, _sieve_prime_cnt@PAGEOFF
    str  w22, [x0]

    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _build_recip_tables - Granlund-Montgomery magic numbers for every
// base prime: pinv (p^-1 mod 2^64), plim ((2^64-1)/p), psq (p^2).
// =====================================================================
_build_recip_tables:
    stp x29, x30, [sp, #-16]!
    mov x29, sp

    adrp x0, _base_prime_cnt@PAGE
    add  x0, x0, _base_prime_cnt@PAGEOFF
    ldr  w19, [x0]                  // count
    adrp x20, _base_primes@PAGE
    add  x20, x20, _base_primes@PAGEOFF
    adrp x21, _pinv_tab@PAGE
    add  x21, x21, _pinv_tab@PAGEOFF
    adrp x22, _plim_tab@PAGE
    add  x22, x22, _plim_tab@PAGEOFF
    adrp x23, _psq_tab@PAGE
    add  x23, x23, _psq_tab@PAGEOFF

    mov  x24, #0                    // i
Lrt_loop:
    cmp  x24, x19
    b.hs Lrt_done

    ldr  w9, [x20, x24, lsl #2]     // p (zero-extended)
    mov  x10, x9                    // x = p (seed)
    mov  x11, #5                    // 5 Newton doublings
Lrt_newton:
    mul  x12, x9, x10               // p*x mod 2^64
    mov  x13, #2
    sub  x13, x13, x12              // 2 - p*x
    mul  x10, x10, x13              // x = x*(2-p*x)
    subs x11, x11, #1
    b.ne Lrt_newton
    str  x10, [x21, x24, lsl #3]    // pinv_tab[i]

    mov  x0, #-1
    udiv x1, x0, x9
    str  x1, [x22, x24, lsl #3]     // plim_tab[i]

    mul  x2, x9, x9
    str  x2, [x23, x24, lsl #3]     // psq_tab[i]

    add  x24, x24, #1
    b    Lrt_loop
Lrt_done:
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _run_wave1 - spawn/join num_threads _SieveThread workers
// =====================================================================
_run_wave1:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!

    adrp x19, _thread_handles@PAGE
    add  x19, x19, _thread_handles@PAGEOFF
    mov  x20, #0
Lw1_spawn:
    add  x0, x19, x20, lsl #3
    mov  x1, #0
    adrp x2, _SieveThread@PAGE
    add  x2, x2, _SieveThread@PAGEOFF
    mov  x3, x20
    bl   _pthread_create
    add  x20, x20, #1
    cmp  x20, #num_threads
    b.lt Lw1_spawn

    mov  x20, #0
Lw1_join:
    add  x0, x19, x20, lsl #3
    ldr  x0, [x0]
    mov  x1, #0
    bl   _pthread_join
    add  x20, x20, #1
    cmp  x20, #num_threads
    b.lt Lw1_join

    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _SieveThread(x0=thread_id) - Wave 1: segmented odds-only prime sieve.
// Persistent regs: x19=seg buf ptr, x20=idx array ptr, x21=local prime
// count, x22=local largest prime, x23=base_prime_cnt, x24=base_primes
// ptr, x26=chunk start, x27=chunk end (exclusive).
// =====================================================================
_SieveThread:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    stp x23, x24, [sp, #-16]!
    stp x25, x26, [sp, #-16]!
    stp x27, x28, [sp, #-16]!

    mov x28, x0                      // thread_id

    adrp x19, _segment_buffers@PAGE
    add  x19, x19, _segment_buffers@PAGEOFF
    mov  x9, #segment_size
    madd x19, x28, x9, x19           // seg buf = base + tid*segment_size

    adrp x20, _next_indices@PAGE
    add  x20, x20, _next_indices@PAGEOFF
    movz x9, #(per_thread_idx_stride & 0xFFFF)
    movk x9, #((per_thread_idx_stride >> 16) & 0xFFFF), lsl #16
    madd x20, x28, x9, x20           // idx arr = base + tid*per_thread_idx_stride

    mov  x21, #0                     // local prime count
    mov  x22, #0                     // local largest prime

    adrp x0, _base_prime_cnt@PAGE
    add  x0, x0, _base_prime_cnt@PAGEOFF
    ldr  w23, [x0]
    adrp x24, _base_primes@PAGE
    add  x24, x24, _base_primes@PAGEOFF

Lsv_get_task:
    adrp x0, _global_task_start@PAGE
    add  x0, x0, _global_task_start@PAGEOFF
    mov  x1, #chunk_size
    ldaddal x1, x26, [x0]             // x26 = old value; [x0] += chunk_size

    LDIMM64 x9, limit
    cmp  x26, x9
    b.hs Lsv_thread_end

    add  x27, x26, #chunk_size
    cmp  x27, x9
    b.ls Lsv_init_indices
    add  x27, x9, #1
Lsv_init_indices:
    mov  x0, #0                      // rbx (prime index)
Lsv_idx_loop:
    cmp  x0, x23
    b.hs Lsv_idx_done
    ldr  w1, [x24, x0, lsl #2]       // p
    mul  x2, x1, x1                   // p*p
    cmp  x2, x26
    b.hs Lsv_idx_store

    udiv x3, x26, x1
    mul  x4, x3, x1
    cmp  x4, x26
    b.hs Lsv_idx_chkodd
    add  x4, x4, x1
Lsv_idx_chkodd:
    tst  x4, #1
    b.ne Lsv_idx_setval
    add  x4, x4, x1
Lsv_idx_setval:
    mov  x2, x4
Lsv_idx_store:
    str  x2, [x20, x0, lsl #3]
    add  x0, x0, #1
    b    Lsv_idx_loop
Lsv_idx_done:

    // memset segment buffer to all-1s (16384 qwords = 131072 bytes)
    mov x0, x19
    mov x1, #16384
    mov x2, #-1
Lsv_fillbuf:
    str x2, [x0], #8
    subs x1, x1, #1
    b.ne Lsv_fillbuf

    cbnz x26, Lsv_sieve
    ldrb w0, [x19]
    and  w0, w0, #0xFE                // n=1 is not prime
    strb w0, [x19]
    add  x21, x21, #1                  // account for prime 2 (not represented)
Lsv_sieve:
    mov  x0, #0                       // rbx
Lsv_prime_loop:
    cmp  x0, x23
    b.hs Lsv_prime_done
    ldr  x9, [x20, x0, lsl #3]        // starting number value
    ldr  w1, [x24, x0, lsl #2]        // p
    cmp  x9, x27
    b.hs Lsv_prime_skip

    lsl  x2, x1, #1                    // step = 2p
Lsv_clear_loop:
    sub  x3, x9, x26
    lsr  x3, x3, #1                    // bit index k
    lsr  x4, x3, #3
    and  x5, x3, #7
    ldrb w6, [x19, x4]
    mov  w7, #1
    lsl  w7, w7, w5
    bic  w6, w6, w7
    strb w6, [x19, x4]
    add  x9, x9, x2
    cmp  x9, x27
    b.lo Lsv_clear_loop
Lsv_prime_skip:
    add x0, x0, #1
    b   Lsv_prime_loop
Lsv_prime_done:

    // ---- NEON popcount over the 131072-byte buffer ----
    movi v31.4s, #0                    // 32-bit accumulator
    mov  x0, x19                        // cursor
    mov  x10, #4                         // outer flush groups
Lsv_pc_outer:
    movi v30.8h, #0                     // 16-bit accumulator (flushed each group)
    mov  x11, #2048                      // 2048*16B = 32768B per group
Lsv_pc_inner:
    ld1  {v0.16b}, [x0], #16
    cnt  v0.16b, v0.16b
    uadalp v30.8h, v0.16b
    subs x11, x11, #1
    b.ne Lsv_pc_inner
    uadalp v31.4s, v30.8h
    subs x10, x10, #1
    b.ne Lsv_pc_outer
    uaddlv d0, v31.4s
    fmov x1, d0
    add  x21, x21, x1

    // ---- find largest prime in this chunk (backward scan) ----
    mov  x9, #131071
    add  x0, x19, x9
    mov  x1, #131072
Lsv_fl_loop:
    ldrb w2, [x0]
    cbnz w2, Lsv_fl_found
    sub  x0, x0, #1
    subs x1, x1, #1
    b.ne Lsv_fl_loop
    b    Lsv_get_task
Lsv_fl_found:
    clz  w3, w2
    mov  w4, #31
    sub  w4, w4, w3                     // bitpos within byte
    sub  x6, x1, #1                      // byte_idx
    lsl  x6, x6, #4
    lsl  x4, x4, #1
    add  x5, x6, x4
    add  x5, x5, x26
    add  x5, x5, #1
    cmp  x5, x27
    b.hs Lsv_fl_retry
    cmp  x5, x22
    b.ls Lsv_get_task
    mov  x22, x5
    b    Lsv_get_task
Lsv_fl_retry:
    sub  x0, x0, #1
    subs x1, x1, #1
    b.ne Lsv_fl_loop
    b    Lsv_get_task

Lsv_thread_end:
    adrp x0, _global_prime_cnt@PAGE
    add  x0, x0, _global_prime_cnt@PAGEOFF
    ldaddal x21, x1, [x0]

    adrp x0, _global_largest@PAGE
    add  x0, x0, _global_largest@PAGEOFF
Lsv_upd_largest:
    ldr   w2, [x0]
    cmp   w22, w2
    b.ls  Lsv_done
    mov   w3, w2
    mov   w4, w22
    casal w3, w4, [x0]
    cmp   w3, w2
    b.ne  Lsv_upd_largest
Lsv_done:
    mov x0, #0
    ldp x27, x28, [sp], #16
    ldp x25, x26, [sp], #16
    ldp x23, x24, [sp], #16
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _run_wave2 - spawn/join num_threads _AmicableThread workers
// =====================================================================
_run_wave2:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!

    adrp x19, _thread_handles@PAGE
    add  x19, x19, _thread_handles@PAGEOFF
    mov  x20, #0
Lw2_spawn:
    add  x0, x19, x20, lsl #3
    mov  x1, #0
    adrp x2, _AmicableThread@PAGE
    add  x2, x2, _AmicableThread@PAGEOFF
    mov  x3, x20
    bl   _pthread_create
    add  x20, x20, #1
    cmp  x20, #num_threads
    b.lt Lw2_spawn

    mov  x20, #0
Lw2_join:
    add  x0, x19, x20, lsl #3
    ldr  x0, [x0]
    mov  x1, #0
    bl   _pthread_join
    add  x20, x20, #1
    cmp  x20, #num_threads
    b.lt Lw2_join

    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _print_wave2_summary(x0=pair_count u64)
// =====================================================================
_print_wave2_summary:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    sub sp, sp, #256

    mov x21, x0
    mov x19, sp
    mov x20, sp

    adrp x1, msg_summary_1@PAGE
    add  x1, x1, msg_summary_1@PAGEOFF
    mov  x0, x20
    mov  x2, #msg_summary_1_len
    bl   _write_str
    mov  x20, x0

    LDIMM64 x0, limit
    mov  x1, x20
    bl   _u64_to_dec
    mov  x20, x0

    adrp x1, msg_summary_2@PAGE
    add  x1, x1, msg_summary_2@PAGEOFF
    mov  x0, x20
    mov  x2, #msg_summary_2_len
    bl   _write_str
    mov  x20, x0

    mov  x0, x21
    mov  x1, x20
    bl   _u64_to_dec
    mov  x20, x0

    adrp x1, msg_summary_3@PAGE
    add  x1, x1, msg_summary_3@PAGEOFF
    mov  x0, x20
    mov  x2, #msg_summary_3_len
    bl   _write_str
    mov  x20, x0

    sub x2, x20, x19
    mov x1, x19
    mov x0, #1
    bl  _write

    add sp, sp, #256
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _FactorizeNumber(x0=N, x1=p_array u64[], x2=e_array u32[], x3=cnt_ptr)
// Trial division (up to sqrt(remaining)); leftover >1 is a final prime.
// At most 9 distinct prime factors are possible for N < 2^32, so the
// 16-slot caller-provided arrays are always large enough.
// =====================================================================
_FactorizeNumber:
    mov x8, x0                  // remaining
    mov x9, #0                    // count

    mov x10, #0                     // exponent of 2
Lfn_div2:
    tst x8, #1
    b.ne Lfn_div2_done
    lsr x8, x8, #1
    add x10, x10, #1
    b Lfn_div2
Lfn_div2_done:
    cbz x10, Lfn_loop3_init
    mov x4, #2
    str x4, [x1, x9, lsl #3]
    str w10, [x2, x9, lsl #2]
    add x9, x9, #1
Lfn_loop3_init:
    mov x11, #3
Lfn_loop3:
    mul x4, x11, x11
    cmp x4, x8
    b.hi Lfn_loop3_done

    udiv x5, x8, x11
    msub x6, x5, x11, x8
    cbnz x6, Lfn_next_p

    mov x10, #0
Lfn_divp:
    udiv x5, x8, x11
    msub x6, x5, x11, x8
    cbnz x6, Lfn_divp_done
    mov x8, x5
    add x10, x10, #1
    b Lfn_divp
Lfn_divp_done:
    str x11, [x1, x9, lsl #3]
    str w10, [x2, x9, lsl #2]
    add x9, x9, #1
Lfn_next_p:
    add x11, x11, #2
    b Lfn_loop3
Lfn_loop3_done:
    cmp x8, #1
    b.ls Lfn_finish
    str x8, [x1, x9, lsl #3]
    mov w4, #1
    str w4, [x2, x9, lsl #2]
    add x9, x9, #1
Lfn_finish:
    str w9, [x3]
    ret

// =====================================================================
// _PrintAmicablePair(x0=N1, x1=N2) - factorizes both, computes the
// Erdos (i,j) classification + regular/irregular ('X') flag, formats
// "[X]i,j\nN1=p^e*...\nN2=p^e*...\n\n" and emits it as a single write(2)
// under print_lock.
//
// Local scratch layout (offset from x21, the sp-relative scratch base):
//   0    p1[16] u64      128  e1[16] u32
//   192  p2[16] u64      320  e2[16] u32
//   384  text buffer (512 bytes)
//   896  cnt_tmp (u32, reused for both factorize calls)
// =====================================================================
_PrintAmicablePair:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    stp x23, x24, [sp, #-16]!
    stp x25, x26, [sp, #-16]!
    stp x27, x28, [sp, #-16]!
    sub sp, sp, #912

    mov x19, x0                 // N1
    mov x20, x1                  // N2
    mov x21, sp                    // scratch base

    mov x0, x19
    add x1, x21, #0
    add x2, x21, #128
    add x3, x21, #896
    bl _FactorizeNumber
    ldr w23, [x21, #896]           // c1

    mov x0, x20
    add x1, x21, #192
    add x2, x21, #320
    add x3, x21, #896
    bl _FactorizeNumber
    ldr w24, [x21, #896]           // c2

    // exponent of prime 2 in each factorization (0 if absent)
    mov x27, #0
    cbz w23, Lpp_chk2_2
    ldr x0, [x21, #0]
    cmp x0, #2
    b.ne Lpp_chk2_2
    ldr w27, [x21, #128]
Lpp_chk2_2:
    mov x28, #0
    cbz w24, Lpp_calc_ij
    ldr x0, [x21, #192]
    cmp x0, #2
    b.ne Lpp_calc_ij
    ldr w28, [x21, #320]

    // Erdos (i,j): merge-walk both (already-sorted) factor lists
Lpp_calc_ij:
    mov x25, #0                     // i
    mov x26, #0                      // j
    mov x0, #0                         // idx1
    mov x1, #0                          // idx2
Lpp_loop_ij:
    cmp w0, w23
    b.lo Lpp_get_pr1
    cmp w1, w24
    b.lo Lpp_get_pr1
    b Lpp_done_ij
Lpp_get_pr1:
    cmp w0, w23
    b.lo Lpp_load_pr1
    mov x2, #-1
    b Lpp_get_pr2
Lpp_load_pr1:
    add x4, x21, #0
    ldr x2, [x4, x0, lsl #3]
Lpp_get_pr2:
    cmp w1, w24
    b.lo Lpp_load_pr2
    mov x3, #-1
    b Lpp_cmp_prs
Lpp_load_pr2:
    add x4, x21, #192
    ldr x3, [x4, x1, lsl #3]
Lpp_cmp_prs:
    cmp x2, x3
    b.eq Lpp_pr_eq
    b.hi Lpp_pr_gt
Lpp_pr_lt:
    add x25, x25, #1
    add x0, x0, #1
    b Lpp_loop_ij
Lpp_pr_gt:
    add x26, x26, #1
    add x1, x1, #1
    b Lpp_loop_ij
Lpp_pr_eq:
    add x4, x21, #128
    ldr w5, [x4, x0, lsl #2]
    add x4, x21, #320
    ldr w6, [x4, x1, lsl #2]
    cmp w5, w6
    b.eq Lpp_adv_both
    b.hi Lpp_inc_i
    add x26, x26, #1
    b Lpp_adv_both
Lpp_inc_i:
    add x25, x25, #1
Lpp_adv_both:
    add x0, x0, #1
    add x1, x1, #1
    b Lpp_loop_ij
Lpp_done_ij:

    add x22, x21, #384              // text cursor

    cmp w27, w28
    b.eq Lpp_fmt_reg
    adrp x1, xflag@PAGE
    add x1, x1, xflag@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0
Lpp_fmt_reg:
    mov x0, x25
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
    adrp x1, comma@PAGE
    add x1, x1, comma@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0
    mov x0, x26
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
    adrp x1, nl@PAGE
    add x1, x1, nl@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    // N1=factors
    mov x0, x19
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
    adrp x1, eqsign@PAGE
    add x1, x1, eqsign@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    mov x25, #0
Lpp_loop_n1:
    cmp w25, w23
    b.hs Lpp_done_n1
    cbz x25, Lpp_no_star1
    adrp x1, star@PAGE
    add x1, x1, star@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0
Lpp_no_star1:
    add x4, x21, #0
    ldr x0, [x4, x25, lsl #3]
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0

    add x4, x21, #128
    ldr w0, [x4, x25, lsl #2]
    cmp w0, #1
    b.ls Lpp_next1

    adrp x1, caret@PAGE
    add x1, x1, caret@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    add x4, x21, #128
    ldr w0, [x4, x25, lsl #2]
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
Lpp_next1:
    add x25, x25, #1
    b Lpp_loop_n1
Lpp_done_n1:

    adrp x1, nl@PAGE
    add x1, x1, nl@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    // N2=factors
    mov x0, x20
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
    adrp x1, eqsign@PAGE
    add x1, x1, eqsign@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    mov x26, #0
Lpp_loop_n2:
    cmp w26, w24
    b.hs Lpp_done_n2
    cbz x26, Lpp_no_star2
    adrp x1, star@PAGE
    add x1, x1, star@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0
Lpp_no_star2:
    add x4, x21, #192
    ldr x0, [x4, x26, lsl #3]
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0

    add x4, x21, #320
    ldr w0, [x4, x26, lsl #2]
    cmp w0, #1
    b.ls Lpp_next2

    adrp x1, caret@PAGE
    add x1, x1, caret@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    add x4, x21, #320
    ldr w0, [x4, x26, lsl #2]
    mov x1, x22
    bl _u64_to_dec
    mov x22, x0
Lpp_next2:
    add x26, x26, #1
    b Lpp_loop_n2
Lpp_done_n2:

    adrp x1, nl@PAGE
    add x1, x1, nl@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0
    adrp x1, nl@PAGE
    add x1, x1, nl@PAGEOFF
    mov x0, x22
    mov x2, #1
    bl _write_str
    mov x22, x0

    // thread-safe emit under print_lock (test-and-set spinlock)
    adrp x9, _print_lock@PAGE
    add x9, x9, _print_lock@PAGEOFF
    mov w10, #1
Lpp_spin:
    ldsetal w10, w11, [x9]
    tst w11, #1
    b.ne Lpp_spin

    add x1, x21, #384
    sub x2, x22, x1
    mov x0, #1
    bl _write

    adrp x9, _print_lock@PAGE
    add x9, x9, _print_lock@PAGEOFF
    stlr wzr, [x9]

    add sp, sp, #912
    ldp x27, x28, [sp], #16
    ldp x25, x26, [sp], #16
    ldp x23, x24, [sp], #16
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _u64_to_dec(x0=value, x1=dest) -> x0=dest+len (no terminator written)
// =====================================================================
_u64_to_dec:
    stp x29, x30, [sp, #-48]!
    mov x29, sp
    mov x2, x0
    add x3, x29, #16
    mov x4, x3
    mov x5, #10
Lud_loop:
    udiv x6, x2, x5
    msub x7, x6, x5, x2
    add  w7, w7, #'0'
    strb w7, [x3], #1
    mov  x2, x6
    cbnz x2, Lud_loop

    mov x8, x1
Lud_rev:
    sub x3, x3, #1
    ldrb w9, [x3]
    strb w9, [x8], #1
    cmp x3, x4
    b.hi Lud_rev
    mov x0, x8
    ldp x29, x30, [sp], #48
    ret

// =====================================================================
// _write_str(x0=dest, x1=src, x2=len) -> x0=dest+len
// =====================================================================
_write_str:
    mov x3, x0
Lws_loop:
    cbz x2, Lws_done
    ldrb w4, [x1], #1
    strb w4, [x3], #1
    sub x2, x2, #1
    b Lws_loop
Lws_done:
    mov x0, x3
    ret

// =====================================================================
// _print_wave1_summary(x0=prime_count u64, x1=largest u32)
// "Primes found: <n>, Largest: <n>\n" -> single write(2) call
// =====================================================================
_print_wave1_summary:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    sub sp, sp, #256

    mov x21, x0
    mov x22, x1

    mov x19, sp
    mov x20, sp

    adrp x1, msg_primes_1@PAGE
    add  x1, x1, msg_primes_1@PAGEOFF
    mov  x0, x20
    mov  x2, #msg_primes_1_len
    bl   _write_str
    mov  x20, x0

    mov  x0, x21
    mov  x1, x20
    bl   _u64_to_dec
    mov  x20, x0

    adrp x1, msg_primes_2@PAGE
    add  x1, x1, msg_primes_2@PAGEOFF
    mov  x0, x20
    mov  x2, #msg_primes_2_len
    bl   _write_str
    mov  x20, x0

    mov  x0, x22
    mov  x1, x20
    bl   _u64_to_dec
    mov  x20, x0

    adrp x1, nl@PAGE
    add  x1, x1, nl@PAGEOFF
    mov  x0, x20
    mov  x2, #1
    bl   _write_str
    mov  x20, x0

    sub x2, x20, x19
    mov x1, x19
    mov x0, #1
    bl  _write

    add sp, sp, #256
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _AmicableThread(x0=thread_id) - Wave 2: segmented divisor sieve +
// memory-lookup amicable pair scan.
// Persistent regs: x19=chunk_rem ptr, x20=chunk_sigma ptr,
// x21=carried-offset array ptr, x22=local pair count,
// x23=sieve_prime_cnt, x24=base_primes ptr, x25=rsi (chunk start),
// x26=chunk end (exclusive), x27/x28=window start/end (relative to
// chunk start); x27 is reused as the am-scan index after windowing.
// =====================================================================
_AmicableThread:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!
    stp x23, x24, [sp, #-16]!
    stp x25, x26, [sp, #-16]!
    stp x27, x28, [sp, #-16]!

    mov x9, x0                        // thread_id

    adrp x19, _chunk_rem@PAGE
    add  x19, x19, _chunk_rem@PAGEOFF
    mov  x10, #(chunk_size*4)
    madd x19, x9, x10, x19

    adrp x20, _chunk_sigma@PAGE
    add  x20, x20, _chunk_sigma@PAGEOFF
    mov  x10, #(chunk_size*8)
    madd x20, x9, x10, x20

    adrp x21, _next_indices@PAGE
    add  x21, x21, _next_indices@PAGEOFF
    movz x10, #(per_thread_idx_stride & 0xFFFF)
    movk x10, #((per_thread_idx_stride >> 16) & 0xFFFF), lsl #16
    madd x21, x9, x10, x21

    mov  x22, #0

    adrp x0, _sieve_prime_cnt@PAGE
    add  x0, x0, _sieve_prime_cnt@PAGEOFF
    ldr  w23, [x0]
    adrp x24, _base_primes@PAGE
    add  x24, x24, _base_primes@PAGEOFF

Lam_get_task:
    adrp x0, _global_task_start@PAGE
    add  x0, x0, _global_task_start@PAGEOFF
    mov  x1, #chunk_size
    ldaddal x1, x25, [x0]

    LDIMM64 x9, limit
    cmp  x25, x9
    b.hs Lam_thread_end

    add  x26, x25, #chunk_size
    cmp  x26, x9
    b.ls Lam_seg_ok
    add  x26, x9, #1
Lam_seg_ok:

    // 1. First-multiple offset of every sieve prime in this chunk
    mov  x0, #0
Lam_calc_off_loop:
    cmp  x0, x23
    b.hs Lam_calc_off_done
    ldr  w1, [x24, x0, lsl #2]        // p
    udiv w2, w25, w1
    msub w3, w2, w1, w25               // remainder = rsi mod p
    mov  w4, #0
    cbz  w3, Lam_off_chk0
    sub  w4, w1, w3
Lam_off_chk0:
    cbnz x25, Lam_off_store
    mov  w4, w1                          // chunk 0: skip N=0, begin at N=p
Lam_off_store:
    str  x4, [x21, x0, lsl #3]
    add  x0, x0, #1
    b    Lam_calc_off_loop
Lam_calc_off_done:

    // 2. Divisor sieve, one L2-resident window (sub_seg_elems) at a time
    mov x27, #0                        // window_start (relative to chunk)
Lam_window_loop:
    sub x9, x26, x25                     // seg_len
    cmp x27, x9
    b.hs Lam_am_scan
    add x28, x27, #sub_seg_elems

    // 2a. sigma[window] = 1
    add x0, x20, x27, lsl #3
    mov x1, #sub_seg_elems
    mov x2, #1
Lam_init_sigma:
    str x2, [x0], #8
    subs x1, x1, #1
    b.ne Lam_init_sigma

    // 2b. rem[window] = N (NEON, 4 lanes/store)
    add x0, x19, x27, lsl #2
    add w1, w25, w27
    dup v0.4s, w1
    adrp x2, iota4@PAGE
    add x2, x2, iota4@PAGEOFF
    ld1 {v1.4s}, [x2]
    add v0.4s, v0.4s, v1.4s
    movi v2.4s, #4
    mov x3, #(sub_seg_elems/4)
Lam_init_rem:
    str q0, [x0], #16
    add v0.4s, v0.4s, v2.4s
    subs x3, x3, #1
    b.ne Lam_init_rem

    // 2c. strip powers of two (even indices only: even N)
    mov x0, x27
Lam_win2_loop:
    cmp x0, x28
    b.hs Lam_win2_done
    add x1, x19, x0, lsl #2
    ldr w2, [x1]
    cbz w2, Lam_win2_next

    mov w3, #2
    mov w4, #3
    lsr w2, w2, #1
    tst w2, #1
    b.ne Lam_win2_store
Lam_win2_while:
    lsr w2, w2, #1
    lsl w3, w3, #1
    add w4, w4, w3
    tst w2, #1
    b.eq Lam_win2_while
Lam_win2_store:
    str w2, [x1]
    add x5, x20, x0, lsl #3
    ldr x6, [x5]
    mul x6, x6, x4
    str x6, [x5]
Lam_win2_next:
    add x0, x0, #2
    b Lam_win2_loop
Lam_win2_done:

    // 2d. odd primes: exact division by pinv, divisibility test vs plim
    mov x0, #0
Lam_wprime_loop:
    cmp x0, x23
    b.hs Lam_wfinish
    ldr x1, [x21, x0, lsl #3]          // carried offset
    cmp x1, x28
    b.hs Lam_wnext_prime

    ldr w2, [x24, x0, lsl #2]           // p (zero-extends x2)
    adrp x3, _pinv_tab@PAGE
    add  x3, x3, _pinv_tab@PAGEOFF
    ldr  x3, [x3, x0, lsl #3]
    adrp x4, _plim_tab@PAGE
    add  x4, x4, _plim_tab@PAGEOFF
    ldr  x4, [x4, x0, lsl #3]
Lam_welem_loop:
    add  x5, x19, x1, lsl #2
    ldr  w6, [x5]
    cbz  w6, Lam_wadvance
    mul  x7, x6, x3                       // rem/p (guaranteed exact)
    mov  x8, x7
    add  x11, x2, #1                        // 1+p
    mov  x10, x2                             // p^1
Lam_wstrip:
    mul  x7, x7, x3
    cmp  x7, x4
    b.hi Lam_wstrip_done
    mov  x8, x7
    mul  x10, x10, x2
    add  x11, x11, x10
    b    Lam_wstrip
Lam_wstrip_done:
    str  w8, [x5]
    add  x12, x20, x1, lsl #3
    ldr  x13, [x12]
    mul  x13, x13, x11
    str  x13, [x12]
Lam_wadvance:
    add  x1, x1, x2
    cmp  x1, x28
    b.lo Lam_welem_loop
    str  x1, [x21, x0, lsl #3]
Lam_wnext_prime:
    add x0, x0, #1
    b Lam_wprime_loop

    // 2e. finalize S(N) for the window while it is still cached
Lam_wfinish:
    mov x0, x27
Lam_wfin_loop:
    cmp x0, x28
    b.hs Lam_wfin_done
    add w1, w25, w0
    cmp w1, #1
    b.ls Lam_wfin_next

    add x2, x19, x0, lsl #2
    ldr w3, [x2]
    add x4, x20, x0, lsl #3
    ldr x5, [x4]

    cmp w3, #1
    b.ls Lam_wfin_norem
    add x6, x3, #1
    mul x5, x5, x6
Lam_wfin_norem:
    sub x5, x5, x1
    str x5, [x4]
Lam_wfin_next:
    add x0, x0, #1
    b Lam_wfin_loop
Lam_wfin_done:
    mov x27, x28
    b Lam_window_loop

    // 3. Instant memory-lookup amicable scan over the whole chunk
Lam_am_scan:
    mov x27, #0
Lam_am_scan_loop:
    sub x9, x26, x25
    cmp x27, x9
    b.hs Lam_get_task

    add w0, w25, w27                     // N
    cmp w0, #1
    b.ls Lam_am_skip

    add x1, x20, x27, lsl #3
    ldr x2, [x1]                           // M = S(N)
    cmp x2, x0
    b.ls Lam_am_skip

    sub x3, x2, x25                          // M - rsi
    cmp x3, x9
    b.hs Lam_am_compute_ssn

    add x4, x20, x3, lsl #3
    ldr x6, [x4]                               // S(M) via direct lookup
    b Lam_am_check
Lam_am_compute_ssn:
    add x1, x2, x0                               // bound = N+M
    str x2, [sp, #-16]!
    mov x0, x2
    bl _SumProperDivisors
    mov x6, x0
    ldr x2, [sp], #16
    add w0, w25, w27                               // recompute N (x0 clobbered by call)
Lam_am_check:
    cmp x6, x0
    b.ne Lam_am_skip

    add x22, x22, #1
    mov x1, x2
    bl _PrintAmicablePair
Lam_am_skip:
    add x27, x27, #1
    b Lam_am_scan_loop

Lam_thread_end:
    adrp x0, _global_pair_cnt@PAGE
    add  x0, x0, _global_pair_cnt@PAGEOFF
    ldaddal x22, x1, [x0]

    mov x0, #0
    ldp x27, x28, [sp], #16
    ldp x25, x26, [sp], #16
    ldp x23, x24, [sp], #16
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// =====================================================================
// _SumProperDivisors(x0=N, x1=abort_bound) -> x0=S(N), or 0 if the
// running sigma product exceeds the bound (no amicable match possible).
// Overflow checks present in the Linux x86 version are omitted: for
// N < 2^32, sigma(N) is bounded well within 2^64 (abundancy index is
// tiny), so the overflow paths are provably unreachable here.
// =====================================================================
_SumProperDivisors:
    stp x29, x30, [sp, #-16]!
    mov x29, sp
    stp x19, x20, [sp, #-16]!
    stp x21, x22, [sp, #-16]!

    mov x19, x1                       // bound
    mov x20, x0                        // remaining cofactor
    mov x22, x0                         // original N
    mov x21, #1                          // sigma accumulator

    tst x20, #1
    b.ne Lsd_skip2

    mov x8, #2
    mov x9, #3
Lsd_div2:
    lsr x20, x20, #1
    tst x20, #1
    b.ne Lsd_div2_done
    lsl x8, x8, #1
    add x9, x9, x8
    b Lsd_div2
Lsd_div2_done:
    mov x21, x9
    cmp x21, x19
    b.hi Lsd_no_match
    // sigma(M) = x21 * sigma(cofactor), so a match needs x21 | bound
    // and bound/x21 (the required sigma(cofactor)) >= cofactor+1.
    udiv x2, x19, x21
    msub x3, x2, x21, x19
    cbnz x3, Lsd_no_match
    cmp  x20, #1
    b.ls Lsd_skip2
    add  x4, x20, #1
    cmp  x2, x4
    b.lo Lsd_no_match

Lsd_skip2:
    mov x0, #0
    adrp x10, _base_prime_cnt@PAGE
    add  x10, x10, _base_prime_cnt@PAGEOFF
    ldr  w10, [x10]
    adrp x11, _base_primes@PAGE
    add  x11, x11, _base_primes@PAGEOFF
    adrp x12, _pinv_tab@PAGE
    add  x12, x12, _pinv_tab@PAGEOFF
    adrp x13, _plim_tab@PAGE
    add  x13, x13, _plim_tab@PAGEOFF
    adrp x14, _psq_tab@PAGE
    add  x14, x14, _psq_tab@PAGEOFF
    mov  x15, #128                  // prime index at which to Miller-Rabin
                                    // the remaining cofactor (see Lsd_mr_check)
Lsd_floop:
    cmp x0, x10
    b.hs Lsd_fend
    cmp x0, x15
    b.eq Lsd_mr_check
Lsd_floop_body:
    ldr x2, [x14, x0, lsl #3]
    cmp x20, x2
    b.lo Lsd_fend

    mov x3, x20
    ldr x4, [x12, x0, lsl #3]
    mul x3, x3, x4
    ldr x5, [x13, x0, lsl #3]
    cmp x3, x5
    b.hi Lsd_np

    ldr w6, [x11, x0, lsl #2]
    mov x7, #1
    mov x8, #1
Lsd_strip:
    mov x20, x3
    mul x8, x8, x6
    add x7, x7, x8
    mul x3, x20, x4
    cmp x3, x5
    b.ls Lsd_strip

    mul x21, x21, x7
    cmp x21, x19
    b.hi Lsd_no_match
    // same divisibility + lower-bound abort as after the 2-strip
    udiv x2, x19, x21
    msub x3, x2, x21, x19
    cbnz x3, Lsd_no_match
    cmp  x20, #1
    b.ls Lsd_np
    add  x4, x20, #1
    cmp  x2, x4
    b.lo Lsd_no_match
    cmp  x20, #0x100, lsl #12       // cofactor still > 2^20: re-arm the
    b.ls Lsd_np                     // Miller-Rabin trigger 128 primes ahead
    add  x15, x0, #128
Lsd_np:
    add x0, x0, #1
    b Lsd_floop
Lsd_fend:
    cmp x20, #1
    b.ls Lsd_done
    add x2, x20, #1
    mul x21, x21, x2
Lsd_done:
    sub x0, x21, x22
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret
Lsd_no_match:
    mov x0, #0
    ldp x21, x22, [sp], #16
    ldp x19, x20, [sp], #16
    ldp x29, x30, [sp], #16
    ret

// Reached (once armed) after 128 consecutive base primes failed to divide
// the cofactor: it then very likely has no small factor at all, so settle
// it with one Miller-Rabin test instead of scanning on to sqrt(cofactor).
Lsd_mr_check:
    mov  x15, #-1                   // disarm (fires at most once per arming)
    cmp  x20, #0x100, lsl #12       // small cofactor: rest of scan is cheap
    b.ls Lsd_floop_body
    stp  x0, x10, [sp, #-48]!       // save scan state (caller-saved regs)
    stp  x11, x12, [sp, #16]
    stp  x13, x14, [sp, #32]
    mov  x0, x20
    bl   _mr_is_prime
    cbnz x0, Lsd_mr_prime
    ldp  x13, x14, [sp, #32]
    ldp  x11, x12, [sp, #16]
    ldp  x0, x10, [sp], #48
    mov  x15, #-1                   // _mr_is_prime clobbered x15: re-disarm
    b    Lsd_floop_body
Lsd_mr_prime:
    add  sp, sp, #48                // cofactor is prime: Lsd_fend finishes
    b    Lsd_fend                   // with sigma *= (cofactor+1) as usual

// =====================================================================
// _mr_is_prime(x0 = n; n odd, 2^20 < n < 2^35) -> x0 = 1 if prime, else 0.
// Deterministic Miller-Rabin with bases {2,3,5,7,11}, exact for all
// n < 2,152,302,898,747 (Jaeschke) - far above any cofactor here, since
// sigma(N) < 4N for N < 2^32 bounds the argument below 2^35.
// Arithmetic is Montgomery multiplication mod n (R = 2^64), reusing the
// Newton-iteration inverse trick from _build_recip_tables.
// Leaf function; clobbers x2-x15 only.
// =====================================================================
_mr_is_prime:
    mov  x2, x0                     // n
    mov  x3, x2                     // Newton for n^-1 mod 2^64 (5 doublings)
    mov  x12, #5
Lmr_newton:
    mul  x13, x2, x3
    mov  x14, #2
    sub  x13, x14, x13
    mul  x3, x3, x13
    subs x12, x12, #1
    b.ne Lmr_newton
    neg  x3, x3                     // -(n^-1) mod 2^64, as REDC wants
    mov  x13, #-1                   // x4 = R mod n (= mont(1))
    udiv x14, x13, x2
    msub x4, x14, x2, x13
    add  x4, x4, #1                 // n odd, so never wraps to 0
    sub  x5, x2, x4                 // x5 = mont(n-1)
    sub  x6, x2, #1                 // n-1 = d * 2^s
    rbit x7, x6
    clz  x7, x7                     // s
    lsr  x6, x6, x7                 // d
    movz x8, #0x7532                // bases 2,3,5,7,11, one nibble each,
    movk x8, #0xB, lsl #16          // consumed low nibble first
Lmr_base_loop:
    and  x9, x8, #15                // a
    lsr  x8, x8, #4
    mul  x9, x9, x4                 // mont(a) = a*R mod n (a*R < 2^40)
    udiv x13, x9, x2
    msub x9, x13, x2, x9
    mov  x10, x4                    // x = mont(1)
    mov  x11, x6                    // remaining bits of d
Lmr_pow_loop:
    tbz  x11, #0, Lmr_pow_sq
    mul   x12, x10, x9              // x = REDC(x * base)
    umulh x13, x10, x9
    mul   x14, x12, x3
    umulh x14, x14, x2
    cmp   x12, #1                   // carry = (lo != 0)
    adc   x10, x13, x14
    subs  x12, x10, x2
    csel  x10, x12, x10, hs
Lmr_pow_sq:
    lsr  x11, x11, #1
    cbz  x11, Lmr_pow_done
    mul   x12, x9, x9               // base = REDC(base^2)
    umulh x13, x9, x9
    mul   x14, x12, x3
    umulh x14, x14, x2
    cmp   x12, #1
    adc   x9, x13, x14
    subs  x12, x9, x2
    csel  x9, x12, x9, hs
    b    Lmr_pow_loop
Lmr_pow_done:
    cmp  x10, x4                    // a^d == 1 ?
    b.eq Lmr_base_pass
    cmp  x10, x5                    // a^d == n-1 ?
    b.eq Lmr_base_pass
    sub  x12, x7, #1                // else square s-1 more times,
Lmr_sq_loop:                        // looking for n-1
    cbz  x12, Lmr_composite
    mul   x13, x10, x10
    umulh x14, x10, x10
    mul   x15, x13, x3
    umulh x15, x15, x2
    cmp   x13, #1
    adc   x10, x14, x15
    subs  x13, x10, x2
    csel  x10, x13, x10, hs
    cmp  x10, x5
    b.eq Lmr_base_pass
    sub  x12, x12, #1
    b    Lmr_sq_loop
Lmr_composite:
    mov  x0, #0
    ret
Lmr_base_pass:
    cbnz x8, Lmr_base_loop
    mov  x0, #1
    ret
