; -------------------------------------------------------------------
; amica_fast_lnx_40.asm - Multi-threaded Segmented Sieve for Linux x64
; Searches the full 2^40 range (0 .. 1,099,511,627,775)
; Optimized for 16-thread Intel Core i5-14400F (AVX2)
; Integrates Odds-Only Prime Sieve & Segmented Divisor Sieve
; Added: Real-time Prime Factorization & Erdos Type Classification
; Perf: DIV-free factor stripping via magic-number reciprocals
;       (Granlund-Montgomery), L2-resident divisor-sieve windows with
;       carried per-prime offsets, THP huge pages, AVX2 window init
;
; Differences vs. the 2^32 build (amica_fast_lnx.asm):
;   * 64-bit residue array: N no longer fits in a dword, so chunk_rem
;     holds qwords (256 MB) next to chunk_sigma (256 MB).
;   * Divisor-sieve prime set grows to pi(2^20) = 82,025 primes; the
;     S(N) fallback needs primes to 2^22 (pi = 295,947) because the
;     partner M can reach ~5.7e12 and must still be trial-divided.
;   * Wave 1 now grabs 32 M-number tasks and walks them as 16 carried
;     sub-segments, so the per-prime start-offset DIV is paid once per
;     32 M numbers instead of once per 2 M - at 2^40 that loop would
;     otherwise cost more than the sieve itself.
;   * Wave 2 start offsets use mulhi-by-reciprocal instead of a 64-bit
;     DIV (~35 cycles -> ~5), and are kept as dwords (chunk-relative)
;     to halve the per-window offset-array traffic.
;   * Miller-Rabin uses 7 bases {2,3,5,7,11,13,17}, deterministic below
;     3.4e14; the 5-base set of the 32-bit build is only exact to
;     2.15e12, under the largest cofactor reachable here.
;   * FactorizeNumber walks base_primes instead of every odd number.
;   * Optional "-n" argument skips Wave 1 (prime statistics), which on
;     its own takes about a minute at this range.
;
; Memory: ~540 MB of .bss (demand-paged, THP-advised).
; -------------------------------------------------------------------
default rel

extern printf
extern sprintf
extern fflush
extern pthread_create
extern pthread_join
extern exit

section .data
    limit          equ 1099511627776   ; 2^40, exclusive upper bound
    sqrt_limit     equ 4194304         ; 2^22: base primes for the S(N) fallback
                                       ; (covers cofactors up to 1.7e13)
    sieve_bound    equ 1048576         ; 2^20 = sqrt(limit): primes both sieves use
    segment_size   equ 131072          ; 128 KB bit buffer (fits in L2)
    segment_span   equ 2097152         ; numbers covered by one bit buffer (odds only)
    sieve_task     equ 33554432        ; 2^25 numbers per Wave-1 task = 16 sub-segments
    chunk_elems    equ 2097152         ; 2^21 numbers per Wave-2 task
    num_threads    equ 16              ; Matches Core i5-14400F logical cores
    sub_seg_elems  equ 32768           ; Divisor-sieve window: 256 KB rem + 256 KB
                                       ; sigma stays L2-resident; must divide chunk
    max_base_pi    equ 300000          ; > pi(2^22) = 295,947
    max_sieve_pi   equ 84000           ; > pi(2^20) =  82,025

    ; Output Formatters
    fmt            db "Primes found: %llu, Largest: %llu", 10, 0
    skip_fmt       db "Prime sieve skipped (-n)", 10, 0
    overflow_fmt   db "Lower int to start in 128-bit range: %llu", 10, 0
    summary_fmt    db "Checked pairs in the range 2..1099511627775; found: %lld amicable pairs", 10, 0

    ; Pair Classification Formatters
    fmt_reg        db "%u,%u", 10, 0
    fmt_irreg      db "X%u,%u", 10, 0
    fmt_N_eq       db "%llu=", 0
    fmt_p          db "%llu", 0
    fmt_p_e        db "%llu^%u", 0
    fmt_str        db "%s", 0

    align 32
    popcount_lut   db 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, \
                      0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
    low_mask       db 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, \
                      0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, \
                      0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, \
                      0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F
    clear_masks    db 0xFE, 0xFD, 0xFB, 0xF7, 0xEF, 0xDF, 0xBF, 0x7F

    align 32
    iota4          dq 0, 1, 2, 3
    const4         dq 4

    print_lock     dd 0

section .bss
    alignb 32
    base_primes       resd max_base_pi
    base_prime_cnt    resd 1
    sieve_prime_cnt   resd 1          ; count of primes < 2^20 (divisor sieve bound)

    ; Magic-number division tables (Granlund-Montgomery), one entry per base prime:
    ; pinv_tab[i] = base_primes[i]^-1 mod 2^64,  plim_tab[i] = (2^64-1)/base_primes[i]
    ; For any x < 2^64:  p | x  <=>  x*pinv (mod 2^64) <= plim,
    ; and when divisible that product is exactly x/p. Replaces DIV with IMUL.
    ; Every base prime is odd, so plim_tab[i] is also floor(2^64 / p), which is
    ; exactly the multiplier the mulhi-based remainder in Wave 2 needs.
    alignb 8
    pinv_tab          resq max_base_pi
    plim_tab          resq max_base_pi
    psq_tab           resq max_base_pi    ; p^2: replaces MUL in the trial-division probe

    ; Thread management
    thread_handles    resq num_threads

    ; Statically allocated memory
    alignb 32
    segment_buffers   resb num_threads * segment_size
    alignb 32
    next_indices      resq num_threads * max_sieve_pi   ; Wave 1: absolute next multiple
    alignb 32
    win_offsets       resd num_threads * max_sieve_pi   ; Wave 2: chunk-relative offset
    alignb 32
    base_sieve_buffer resb sqrt_limit / 16

    ; Segmented Divisor Sieve Buffers (16 MB rem + 16 MB sigma per thread)
    alignb 32
    chunk_rem         resq num_threads * chunk_elems
    alignb 32
    chunk_sigma       resq num_threads * chunk_elems

    ; Atomic Dynamic Work Queue & Synchronization
    alignb 8
    global_task_start resq 1
    global_prime_cnt  resq 1
    global_pair_cnt   resq 1
    global_largest    resq 1
    overflow_reported resq 1

section .text
    global main

main:
    push rbp
    mov rbp, rsp
    sub rsp, 80             ; Keep stack 16-byte aligned

    mov [rbp-8], rdi        ; argc
    mov [rbp-16], rsi       ; argv

    mov dword [rel print_lock], 0

    ; Ask the kernel for transparent huge pages on the big sieve buffers
    ; madvise(addr & ~4095, len, MADV_HUGEPAGE) - advisory, errors ignored
    lea rdi, [rel chunk_rem]
    and rdi, -4096
    mov esi, num_threads * chunk_elems * 8 + 4096
    mov edx, 14             ; MADV_HUGEPAGE
    mov eax, 28             ; sys_madvise
    syscall
    lea rdi, [rel chunk_sigma]
    and rdi, -4096
    mov esi, num_threads * chunk_elems * 8 + 4096
    mov edx, 14
    mov eax, 28
    syscall

    ; 1. Setup base sieve (0 to 4,194,304)
    lea rdi, [rel base_sieve_buffer]
    mov rcx, sqrt_limit / 128
    mov rax, 0xFFFFFFFFFFFFFFFF
    rep stosq

    lea rbx, [rel base_sieve_buffer]
    mov byte [rbx], 0xFE    ; Bit 0 is '1', not prime.

    mov rsi, 3
.base_sieve:
    mov rax, rsi
    shr rax, 1
    lea rbx, [rel base_sieve_buffer]
    bt [rbx], rax
    jnc .next_base

    mov rax, rsi
    mul rsi
    mov r8, rax
.clear_base:
    mov rax, r8
    shr rax, 1
    lea rbx, [rel base_sieve_buffer]
    btr [rbx], rax
    lea rdx, [rsi + rsi]
    add r8, rdx
    cmp r8, sqrt_limit
    jb .clear_base
.next_base:
    add rsi, 2
    cmp rsi, 2049           ; 2047^2 < sqrt_limit <= 2049^2
    jb .base_sieve

    ; Collect base primes
    xor r12, r12
    mov rsi, 3
.collect:
    mov rax, rsi
    shr rax, 1
    lea rbx, [rel base_sieve_buffer]
    bt [rbx], rax
    jnc .not_p
    lea rbx, [rel base_primes]
    mov [rbx + r12*4], esi
    inc r12
.not_p:
    add rsi, 2
    cmp rsi, sqrt_limit
    jb .collect
    mov dword [rel base_prime_cnt], r12d

    ; Count primes below 2^20 (only those are used by the two sieves)
    xor rbx, rbx
    lea rdx, [rel base_primes]
.count_sieve:
    cmp ebx, r12d
    jae .count_done
    cmp dword [rdx + rbx*4], sieve_bound
    jae .count_done
    inc rbx
    jmp .count_sieve
.count_done:
    mov dword [rel sieve_prime_cnt], ebx

    ; Build magic-number reciprocal tables for every base prime (all odd)
    xor rbx, rbx
.build_tabs:
    cmp ebx, r12d
    jae .tabs_done
    lea rdx, [rel base_primes]
    mov ecx, [rdx + rbx*4]  ; p (zero-extended into rcx)
    mov rax, rcx            ; Newton seed: x = p is correct to 3 bits (p*p == 1 mod 8)
    mov r8d, 5              ; 5 doublings: 3 -> 96 correct bits >= 64
.newton:
    mov r9, rcx
    imul r9, rax
    mov r10, 2
    sub r10, r9
    imul rax, r10           ; x = x * (2 - p*x)  (mod 2^64)
    dec r8d
    jnz .newton
    lea rdx, [rel pinv_tab]
    mov [rdx + rbx*8], rax
    xor edx, edx
    mov rax, -1
    div rcx                 ; (2^64-1) / p  == floor(2^64 / p) for odd p
    lea rdx, [rel plim_tab]
    mov [rdx + rbx*8], rax
    mov rax, rcx
    imul rax, rcx           ; p^2
    lea rdx, [rel psq_tab]
    mov [rdx + rbx*8], rax
    inc rbx
    jmp .build_tabs
.tabs_done:

    ; Initialize Globals
    mov qword [rel global_task_start], 0
    mov qword [rel global_prime_cnt], 0
    mov qword [rel global_largest], 0
    mov qword [rel global_pair_cnt], 0
    mov qword [rel overflow_reported], 0

    ; -------------------------------------------------------------------
    ; WAVE 1: Prime Sieve (skipped with "-n")
    ; -------------------------------------------------------------------
    mov rax, [rbp-8]
    cmp rax, 2
    jl .run_sieve
    mov rdx, [rbp-16]
    mov rdx, [rdx + 8]      ; argv[1]
    cmp word [rdx], 0x6E2D  ; "-n"
    jne .run_sieve
    lea rdi, [rel skip_fmt]
    xor eax, eax
    call printf
    jmp .wave2

.run_sieve:
    xor rbx, rbx
.spawn_sieve:
    lea rdi, [rel thread_handles + rbx*8]
    xor rsi, rsi
    lea rdx, [rel SieveThread]
    mov rcx, rbx
    call pthread_create

    inc rbx
    cmp rbx, num_threads
    jb .spawn_sieve

    ; Wait for Sieve Threads
    xor rbx, rbx
.join_sieve:
    mov rdi, [rel thread_handles + rbx*8]
    xor rsi, rsi
    call pthread_join
    inc rbx
    cmp rbx, num_threads
    jb .join_sieve

    ; Print Prime Statistics
    lea rdi, [rel fmt]
    mov rsi, [rel global_prime_cnt]
    mov rdx, [rel global_largest]
    xor eax, eax
    call printf

    ; -------------------------------------------------------------------
    ; WAVE 2: Amicable Pair Search (Segmented Divisor Sieve)
    ; -------------------------------------------------------------------
.wave2:
    mov qword [rel global_task_start], 0

    xor rbx, rbx
.spawn_amica:
    lea rdi, [rel thread_handles + rbx*8]
    xor rsi, rsi
    lea rdx, [rel AmicableThread]
    mov rcx, rbx
    call pthread_create

    inc rbx
    cmp rbx, num_threads
    jb .spawn_amica

    ; Wait for Amicable Threads
    xor rbx, rbx
.join_amica:
    mov rdi, [rel thread_handles + rbx*8]
    xor rsi, rsi
    call pthread_join
    inc rbx
    cmp rbx, num_threads
    jb .join_amica

    ; Print Amicable Pairs Summary
    lea rdi, [rel summary_fmt]
    mov rsi, [rel global_pair_cnt]
    xor eax, eax
    call printf

    xor rdi, rdi
    call exit

; -------------------------------------------------------------------
; SieveThread - Pass 1: Identifies Primes
; A task is 32 M numbers walked as 16 carried sub-segments, so the
; per-prime first-multiple DIV is amortized 16x. At 2^40 with 82 K
; sieve primes that loop would otherwise dominate the whole wave.
; Stack: [rsp+0] task_end, [rsp+8] task_start
; -------------------------------------------------------------------
SieveThread:
    push rbp
    push rbx
    push r12
    push r13
    push r14
    push r15
    sub rsp, 72

    mov rbp, rdi

    mov rax, rbp
    shl rax, 17                 ; tid * segment_size
    lea r14, [rel segment_buffers]
    add r14, rax

    mov rax, rbp
    imul rax, max_sieve_pi * 8
    lea r15, [rel next_indices]
    add r15, rax

    xor r12, r12                ; primes found by this thread
    xor r13, r13                ; largest prime seen by this thread

.get_task:
    mov rax, sieve_task
    lea rcx, [rel global_task_start]
    lock xadd [rcx], rax
    mov rsi, rax

    mov rax, limit
    cmp rsi, rax
    jae .thread_end

    mov r10, rsi
    add r10, sieve_task
    cmp r10, rax
    jbe .init_indices
    mov r10, rax

.init_indices:
    mov [rsp + 0], r10          ; task_end
    mov [rsp + 8], rsi          ; task_start

    xor rbx, rbx
    mov ecx, dword [rel sieve_prime_cnt]
.idx_loop:
    lea rdx, [rel base_primes]
    mov eax, [rdx + rbx*4]
    mov r8, rax
    imul rax, rax               ; p^2 < 2^40
    cmp rax, rsi
    jae .store_idx
    mov rax, rsi
    xor rdx, rdx
    div r8
    imul rax, r8                ; largest multiple of p <= task_start
    cmp rax, rsi
    jae .chk_odd_idx
    add rax, r8
.chk_odd_idx:
    test rax, 1
    jnz .store_idx
    add rax, r8                 ; p is odd, so this flips parity
.store_idx:
    mov [r15 + rbx*8], rax
    inc rbx
    cmp ebx, ecx
    jb .idx_loop

    ; Walk the task one bit-buffer sub-segment at a time
.subseg_loop:
    mov r10, rsi
    add r10, segment_span
    cmp r10, [rsp + 0]
    jbe .have_end
    mov r10, [rsp + 0]
.have_end:

    mov rdi, r14
    mov rcx, segment_size / 8
    mov rax, 0xFFFFFFFFFFFFFFFF
    rep stosq

    test rsi, rsi
    jnz .sieve
    btr qword [r14], 0          ; 1 is not prime
    inc r12                     ; 2 is, and the odds-only buffer omits it

.sieve:
    xor rbx, rbx
    mov ecx, dword [rel sieve_prime_cnt]
    lea r8, [rel clear_masks]
    lea r11, [rel base_primes]
.prime_loop:
    mov r9, [r15 + rbx*8]
    cmp r9, r10
    jae .skip_p
    mov eax, [r11 + rbx*4]
    lea rax, [rax + rax]        ; 2p: stay on odd multiples
.clear:
    mov rdi, r9
    sub rdi, rsi
    shr rdi, 1
    mov rdx, rdi
    shr rdi, 3
    and rdx, 7
    mov dl, [r8 + rdx]
    and [r14 + rdi], dl
    add r9, rax
    cmp r9, r10
    jb .clear
    mov [r15 + rbx*8], r9       ; carry into the next sub-segment
.skip_p:
    inc rbx
    cmp ebx, ecx
    jb .prime_loop

    ; Fast AVX2 Count
    vpxor ymm0, ymm0, ymm0
    vpxor ymm6, ymm6, ymm6
    lea rdx, [rel popcount_lut]
    vmovdqu ymm4, [rdx]
    lea rdx, [rel low_mask]
    vmovdqu ymm5, [rdx]
    mov rdi, r14
    mov rcx, segment_size / 32
.count:
    vmovdqu ymm1, [rdi]
    vpand ymm2, ymm1, ymm5
    vpsrlw ymm3, ymm1, 4
    vpand ymm3, ymm3, ymm5
    vpshufb ymm2, ymm4, ymm2
    vpshufb ymm3, ymm4, ymm3
    vpaddb ymm1, ymm2, ymm3
    vpsadbw ymm1, ymm1, ymm6
    vpaddq ymm0, ymm0, ymm1
    add rdi, 32
    dec rcx
    jnz .count
    vextracti128 xmm1, ymm0, 1
    vpaddq xmm0, xmm0, xmm1
    movq rax, xmm0
    psrldq xmm0, 8
    movq rdx, xmm0
    add rax, rdx
    add r12, rax

    ; Find Largest Prime in this sub-segment
    lea rdi, [r14 + segment_size - 1]
    mov rcx, segment_size
.f_last:
    movzx eax, byte [rdi]
    test eax, eax
    jnz .found_p
    dec rdi
    dec rcx
    jnz .f_last
    jmp .next_subseg
.found_p:
    bsr eax, eax
    mov rdx, rcx
    dec rdx                     ; byte index inside the buffer
    shl rdx, 4                  ; 16 numbers per byte
    lea rax, [rax + rax]
    add rax, rdx
    add rax, rsi
    inc rax                     ; N = start + 16*byte + 2*bit + 1
    cmp rax, r10
    jae .retry_last
    cmp rax, r13
    cmova r13, rax
    jmp .next_subseg
.retry_last:
    dec rdi
    dec rcx
    jnz .f_last

.next_subseg:
    add rsi, segment_span
    cmp rsi, [rsp + 0]
    jb .subseg_loop
    mov rsi, [rsp + 8]
    jmp .get_task

.thread_end:
    lea rcx, [rel global_prime_cnt]
    mov rax, r12
    lock xadd [rcx], rax

.upd_largest:
    mov rax, qword [rel global_largest]
    cmp r13, rax
    jbe .done_end
    mov rdx, r13
    lock cmpxchg qword [rel global_largest], rdx
    jnz .upd_largest

.done_end:
    vzeroupper
    add rsp, 72
    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    pop rbp
    xor rax, rax
    ret

; -------------------------------------------------------------------
; AmicableThread - Pass 2: Segmented Divisor Sieve & Pairs Scan
; Processes each 2 M chunk in L2-resident windows of sub_seg_elems
; elements; all factor stripping is DIV-free via pinv/plim tables.
; Stack: [rsp+40] win_offsets base, [rsp+48] window start,
;        [rsp+56] window end, [rsp+64] seg_len
; -------------------------------------------------------------------
AmicableThread:
    push rbp
    push rbx
    push r12
    push r13
    push r14
    push r15
    sub rsp, 72

    mov rbp, rdi

    mov rax, rbp
    shl rax, 21                 ; tid * chunk_elems
    lea r8, [rax * 8]

    lea r14, [rel chunk_rem]
    add r14, r8

    lea r15, [rel chunk_sigma]
    add r15, r8

    ; Per-thread carried-offset array (chunk-relative, dwords)
    mov rax, rbp
    imul rax, max_sieve_pi * 4
    lea rcx, [rel win_offsets]
    add rcx, rax
    mov [rsp + 40], rcx

    xor r12, r12

.get_task:
    mov rax, chunk_elems
    lea rcx, [rel global_task_start]
    lock xadd [rcx], rax
    mov rsi, rax

    mov rax, limit
    cmp rsi, rax
    jae .thread_end

    mov r10, rsi
    add r10, chunk_elems
    cmp r10, rax
    jbe .set_end
    mov r10, rax
.set_end:
    mov rax, r10
    sub rax, rsi
    mov [rsp + 64], rax         ; seg_len (always chunk_elems here)

    ; 1. First-multiple offset of every sieve prime in this chunk.
    ;    q_hat = mulhi(chunk_start, floor(2^64/p)) is either q or q-1,
    ;    so one conditional subtract fixes the remainder - no 64-bit DIV.
    mov ebx, dword [rel sieve_prime_cnt]
    xor r13d, r13d
    mov rcx, [rsp + 40]
    lea rdi, [rel base_primes]
.calc_off:
    cmp r13d, ebx
    jae .offs_done
    mov r9d, [rdi + r13*4]                          ; p
    mov rax, [rdi + (plim_tab - base_primes) + r13*8]
    mul rsi                                         ; rdx = q_hat
    mov rax, rdx
    imul rax, r9
    mov r8, rsi
    sub r8, rax                                     ; remainder in [0, 2p)
    cmp r8, r9
    jb .rem_ok
    sub r8, r9
.rem_ok:
    mov eax, r9d
    sub rax, r8                                     ; p - rem
    test r8, r8
    cmovz rax, r8                                   ; exact multiple -> offset 0
    test rsi, rsi
    jnz .off_store
    mov eax, r9d                                    ; chunk 0: skip N=0, begin at N=p
.off_store:
    mov [rcx + r13*4], eax
    inc r13d
    jmp .calc_off
.offs_done:

    ; 2. Divisor sieve, one L2-resident window at a time
    mov qword [rsp + 48], 0             ; window start (element offset in chunk)
.window_loop:
    mov r11, [rsp + 48]
    cmp r11, [rsp + 64]
    jae .am_scan
    lea rdx, [r11 + sub_seg_elems]
    mov [rsp + 56], rdx                 ; window end

    ; 2a. Initialize Sigma Array to 1
    lea rdi, [r15 + r11*8]
    mov rcx, sub_seg_elems
    mov rax, 1
    rep stosq

    ; 2b. Initialize Remainder Array to N (AVX2, 4 qword lanes per store)
    lea rdi, [r14 + r11*8]
    mov rax, rsi
    add rax, r11
    vmovq xmm0, rax
    vpbroadcastq ymm0, xmm0
    vpaddq ymm0, ymm0, [rel iota4]
    vpbroadcastq ymm1, qword [rel const4]
    mov rcx, sub_seg_elems / 4
.init_rem:
    vmovdqa [rdi], ymm0
    vpaddq ymm0, ymm0, ymm1
    add rdi, 32
    dec rcx
    jnz .init_rem
    vzeroupper

    ; 2c. Strip powers of two (even N; chunk and window starts are even)
    mov r8, r11
.win2_loop:
    cmp r8, [rsp + 56]
    jae .win2_done
    mov rax, [r14 + r8*8]
    test rax, rax
    jz .win2_next

    mov r10, 2
    mov rdx, 3
    shr rax, 1
    mov rcx, rax

    test rax, 1
    jnz .end_div_2
.div_2_while:
    shr rax, 1
    mov rcx, rax
    shl r10, 1
    add rdx, r10
    test rax, 1
    jz .div_2_while
.end_div_2:
    mov [r14 + r8*8], rcx
    mov rax, [r15 + r8*8]
    imul rax, rdx
    mov [r15 + r8*8], rax
.win2_next:
    add r8, 2
    jmp .win2_loop
.win2_done:

    ; 2d. Odd primes: exact division by pinv, divisibility test vs plim
    lea rdi, [rel base_primes]
    xor r13d, r13d
.wprime_loop:
    cmp r13d, ebx
    jae .wfinish
    mov rcx, [rsp + 40]
    mov r8d, [rcx + r13*4]              ; carried offset of next multiple
    cmp r8, [rsp + 56]
    jae .wnext_prime
    mov r9d, [rdi + r13*4]              ; p
    mov rbp, [rdi + (pinv_tab - base_primes) + r13*8]
    mov rdx, [rdi + (plim_tab - base_primes) + r13*8]
.welem_loop:
    mov rax, [r14 + r8*8]
    test rax, rax
    jz .wadvance
    imul rax, rbp                       ; rem / p (exact: p | rem here)
    mov rcx, rax
    lea r11, [r9 + 1]                   ; 1 + p
    mov r10, r9                         ; p^k
.wstrip:
    imul rax, rbp
    cmp rax, rdx
    ja .wstrip_done                     ; p no longer divides
    mov rcx, rax
    imul r10, r9
    add r11, r10
    jmp .wstrip
.wstrip_done:
    mov [r14 + r8*8], rcx
    mov rax, [r15 + r8*8]
    imul rax, r11                       ; sigma < 2^43, no overflow
    mov [r15 + r8*8], rax
.wadvance:
    add r8, r9
    cmp r8, [rsp + 56]
    jb .welem_loop
    mov rcx, [rsp + 40]
    mov [rcx + r13*4], r8d              ; carry into next window
.wnext_prime:
    inc r13d
    jmp .wprime_loop

    ; 2e. Calculate S(N) for the window while it is still cached
.wfinish:
    mov r8, [rsp + 48]
.wfin_loop:
    cmp r8, [rsp + 56]
    jae .wfin_done
    mov rax, rsi
    add rax, r8                         ; N
    cmp rax, 1
    jbe .wfin_next

    mov rcx, [r14 + r8*8]
    mov r11, [r15 + r8*8]

    cmp rcx, 1
    jbe .wfin_norem
    inc rcx
    imul r11, rcx                       ; sigma *= (rem + 1); rem is prime here
.wfin_norem:
    sub r11, rax
    mov [r15 + r8*8], r11
.wfin_next:
    inc r8
    jmp .wfin_loop
.wfin_done:
    mov rax, [rsp + 56]
    mov [rsp + 48], rax
    jmp .window_loop

    ; 3. Instant Memory Lookup Amicable Scan
.am_scan:
    xor rdi, rdi
.am_scan_loop:
    cmp rdi, [rsp + 64]
    jae .get_task

    mov rax, rsi
    add rax, rdi                        ; N
    cmp rax, 1
    jbe .skip_am

    mov r13, [r15 + rdi*8]              ; M = S(N)
    cmp r13, rax
    jbe .skip_am                        ; only report from the smaller member

    mov r9, r13
    sub r9, rsi
    cmp r9, [rsp + 64]
    jae .compute_ssn                    ; M outside this chunk (also catches M < start)

    mov r10, [r15 + r9*8]
    jmp .check_pair

.compute_ssn:
    lea rdx, [r13 + rax]    ; abort bound: a match needs sigma(M) == N + M
    mov rcx, r13
    call SumProperDivisors
    cmp rax, -1
    je .handle_ovf
    mov r10, rax

.check_pair:
    mov rax, rsi
    add rax, rdi
    cmp r10, rax
    jne .skip_am

    ; PAIR FOUND!
    inc r12
    push rdi                ; Save loop index
    push rsi                ; Save chunk_start

    mov rdi, rax            ; Factorize and Format N1
    mov rsi, r13            ; Factorize and Format N2
    call PrintAmicablePair

    pop rsi                 ; Restore chunk_start
    pop rdi                 ; Restore loop index
    jmp .skip_am

.handle_ovf:
    xor rax, rax
    mov rdx, 1
    lock cmpxchg qword [rel overflow_reported], rdx
    jnz .skip_am

    push rdi
    push rsi
    mov rax, rsi
    add rax, rdi
    lea rdi, [rel overflow_fmt]
    mov rsi, rax
    xor eax, eax
    call printf
    pop rsi
    pop rdi

.skip_am:
    inc rdi
    jmp .am_scan_loop

.thread_end:
    lea rcx, [rel global_pair_cnt]
    mov rax, r12
    lock xadd [rcx], rax

    vzeroupper
    add rsp, 72
    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    pop rbp
    xor rax, rax
    ret

; -------------------------------------------------------------------
; PrintAmicablePair - Factorizes and formats pairs with Erdos (i,j)
; -------------------------------------------------------------------
PrintAmicablePair:
    push rbp
    mov rbp, rsp
    sub rsp, 2048       ; Space for local arrays and format buffer
    push rbx
    push r12
    push r13
    push r14
    push r15
    sub rsp, 8          ; Ensure 16-byte stack alignment prior to calls

    ; RDI = N1, RSI = N2
    mov [rbp-8], rdi    ; Save N1
    mov [rbp-16], rsi   ; Save N2

    ; Factorize N1
    lea rsi, [rbp-296]  ; p1 array
    lea rdx, [rbp-424]  ; e1 array
    lea rcx, [rbp-20]   ; c1 size_ptr
    call FactorizeNumber

    ; Factorize N2
    mov rdi, [rbp-16]
    lea rsi, [rbp-680]  ; p2 array
    lea rdx, [rbp-808]  ; e2 array
    lea rcx, [rbp-24]   ; c2 size_ptr
    call FactorizeNumber

    mov dword [rbp-28], 0   ; i
    mov dword [rbp-32], 0   ; j
    mov dword [rbp-36], 0   ; exp2_1
    mov dword [rbp-40], 0   ; exp2_2

    ; Check if power of 2 differs (determines Regular vs Irregular 'X')
    mov eax, [rbp-20]
    test eax, eax
    jz .chk_exp2_2
    cmp qword [rbp-296], 2
    jne .chk_exp2_2
    mov eax, [rbp-424]
    mov [rbp-36], eax
.chk_exp2_2:
    mov eax, [rbp-24]
    test eax, eax
    jz .calc_ij
    cmp qword [rbp-680], 2
    jne .calc_ij
    mov eax, [rbp-808]
    mov [rbp-40], eax

.calc_ij:
    xor r8d, r8d    ; idx1 = 0
    xor r9d, r9d    ; idx2 = 0

.loop_ij:
    mov eax, [rbp-20]
    mov ebx, [rbp-24]
    cmp r8d, eax
    jb .get_pr1
    cmp r9d, ebx
    jb .get_pr1
    jmp .done_ij

.get_pr1:
    cmp r8d, eax
    jb .load_pr1
    mov r10, 0xFFFFFFFFFFFFFFFF ; pr1 = INF
    jmp .get_pr2
.load_pr1:
    mov r10, [rbp-296 + r8*8]

.get_pr2:
    cmp r9d, ebx
    jb .load_pr2
    mov r11, 0xFFFFFFFFFFFFFFFF ; pr2 = INF
    jmp .cmp_prs
.load_pr2:
    mov r11, [rbp-680 + r9*8]

.cmp_prs:
    cmp r10, r11
    je .pr_eq
    ja .pr_gt

.pr_lt:
    inc dword [rbp-28]  ; i++
    inc r8d             ; idx1++
    jmp .loop_ij

.pr_gt:
    inc dword [rbp-32]  ; j++
    inc r9d             ; idx2++
    jmp .loop_ij

.pr_eq:
    mov ecx, [rbp-424 + r8*4]
    mov edx, [rbp-808 + r9*4]
    cmp ecx, edx
    je .adv_both
    ja .inc_i
    inc dword [rbp-32]  ; j++
    jmp .adv_both
.inc_i:
    inc dword [rbp-28]  ; i++
.adv_both:
    inc r8d
    inc r9d
    jmp .loop_ij

.done_ij:
    ; Format into [rbp-2040]
    lea rdi, [rbp-2040]
    mov eax, [rbp-36]
    cmp eax, [rbp-40]
    je .fmt_reg
    lea rsi, [rel fmt_irreg]
    jmp .do_fmt1
.fmt_reg:
    lea rsi, [rel fmt_reg]
.do_fmt1:
    mov edx, [rbp-28]
    mov ecx, [rbp-32]
    xor eax, eax
    call sprintf

    lea r12, [rbp-2040]
    movsxd rax, eax
    add r12, rax

    ; Output string N1=
    mov rdi, r12
    lea rsi, [rel fmt_N_eq]
    mov rdx, [rbp-8]
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax

    xor r13d, r13d           ; k=0
.loop_N1:
    cmp r13d, [rbp-20]
    jge .done_N1

    test r13d, r13d
    jz .no_star1
    mov byte [r12], '*'
    inc r12
    mov byte [r12], 0
.no_star1:
    mov eax, [rbp-424 + r13*4]
    cmp eax, 1
    ja .pow1

    mov rdi, r12
    lea rsi, [rel fmt_p]
    mov rdx, [rbp-296 + r13*8]
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax
    jmp .next1
.pow1:
    mov rdi, r12
    lea rsi, [rel fmt_p_e]
    mov rdx, [rbp-296 + r13*8]
    mov ecx, eax
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax
.next1:
    inc r13d
    jmp .loop_N1
.done_N1:

    mov byte [r12], 10      ; '\n'
    inc r12
    mov byte [r12], 0

    ; Output string N2=
    mov rdi, r12
    lea rsi, [rel fmt_N_eq]
    mov rdx, [rbp-16]
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax

    xor r13d, r13d
.loop_N2:
    cmp r13d, [rbp-24]
    jge .done_N2

    test r13d, r13d
    jz .no_star2
    mov byte [r12], '*'
    inc r12
    mov byte [r12], 0
.no_star2:
    mov eax, [rbp-808 + r13*4]
    cmp eax, 1
    ja .pow2

    mov rdi, r12
    lea rsi, [rel fmt_p]
    mov rdx, [rbp-680 + r13*8]
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax
    jmp .next2
.pow2:
    mov rdi, r12
    lea rsi, [rel fmt_p_e]
    mov rdx, [rbp-680 + r13*8]
    mov ecx, eax
    xor eax, eax
    call sprintf
    movsxd rax, eax
    add r12, rax
.next2:
    inc r13d
    jmp .loop_N2
.done_N2:

    ; Append double newline matching C++ format
    mov byte [r12], 10
    inc r12
    mov byte [r12], 10
    inc r12
    mov byte [r12], 0

    ; Thread-safe Terminal Print
.spin_print:
    lock bts dword [rel print_lock], 0
    jc .spin_print

    lea rdi, [rel fmt_str]
    lea rsi, [rbp-2040]
    xor eax, eax
    call printf

    xor rdi, rdi            ; flush all streams
    call fflush

    mov dword [rel print_lock], 0

    add rsp, 8
    pop r15
    pop r14
    pop r13
    pop r12
    pop rbx
    mov rsp, rbp
    pop rbp
    ret

; -------------------------------------------------------------------
; FactorizeNumber - Generates prime factors and powers
; RDI=N, RSI=p_array, RDX=e_array, RCX=count_ptr
; Walks base_primes (up to 2^22) rather than every odd number - the
; arguments here reach ~5.7e12, so that is 8x fewer trial divisions
; and the table fully covers sqrt(N).
; -------------------------------------------------------------------
FactorizeNumber:
    push rbx
    push r12
    push r13
    push r14

    mov r8, rdi     ; N
    mov r12, rdx    ; e_array
    mov r13, rcx    ; c_ptr
    xor r9, r9      ; count = 0

    xor r10, r10    ; exp = 0
.div2:
    test r8, 1
    jnz .done2
    inc r10
    shr r8, 1
    jmp .div2
.done2:
    test r10, r10
    jz .loop3_init
    mov qword [rsi + r9*8], 2
    mov dword [r12 + r9*4], r10d
    inc r9

.loop3_init:
    xor r14, r14                ; index into base_primes
    lea rbx, [rel base_primes]
    mov ecx, dword [rel base_prime_cnt]
.loop3:
    cmp r14d, ecx
    jae .done3
    mov r11d, [rbx + r14*4]     ; p
    mov rax, r11
    imul rax, r11
    cmp rax, r8
    ja .done3                   ; p^2 > remaining: what is left is prime

    mov rax, r8
    xor rdx, rdx
    div r11
    test rdx, rdx
    jnz .next_p

    xor r10, r10    ; exp = 0
.div_p:
    mov r8, rax     ; N = N / p
    inc r10
    mov rax, r8
    xor rdx, rdx
    div r11
    test rdx, rdx
    jz .div_p

    mov qword [rsi + r9*8], r11
    mov dword [r12 + r9*4], r10d
    inc r9

.next_p:
    inc r14
    jmp .loop3

.done3:
    cmp r8, 1
    jbe .finish
    mov qword [rsi + r9*8], r8
    mov dword [r12 + r9*4], 1
    inc r9
.finish:
    mov dword [r13], r9d

    pop r14
    pop r13
    pop r12
    pop rbx
    ret

; -------------------------------------------------------------------
; SumProperDivisors - Fallback S(N), DIV-free via pinv/plim tables
; RCX = N (preserved), RDX = abort bound (the running sigma product
; can only grow, so once it exceeds the bound no match is possible).
; Returns RAX = S(N), 0 on early abort, or -1 on 64-bit overflow.
; -------------------------------------------------------------------
SumProperDivisors:
    push rbp
    push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    mov rbp, rdx            ; bound = N + M for the amicable check
    mov r12, rcx
    mov r13, 1

    test r12, 1
    jnz .skip_2
    mov r8, 1
    mov r9, 1
.div_2:
    shr r12, 1
    shl r9, 1
    jc .ovf
    add r8, r9
    jc .ovf
    test r12, 1
    jz .div_2
    mov r13, r8
    cmp r13, rbp
    ja .no_match
    ; sigma(M) = r13 * sigma(cofactor), so a match needs r13 | bound
    ; and bound/r13 (the required sigma(cofactor)) >= cofactor+1.
    mov rax, rbp
    xor edx, edx
    div r13
    test rdx, rdx
    jnz .no_match
    cmp r12, 1
    jbe .skip_2
    lea rdx, [r12 + 1]
    cmp rax, rdx
    jb .no_match

.skip_2:
    xor rsi, rsi
    lea r14, [rel base_primes]
    mov edi, dword [rel base_prime_cnt]
    mov r11, 128            ; prime index at which to Miller-Rabin the
                            ; remaining cofactor (see .mr_check)

.f_loop:
    cmp esi, edi
    jae .ovf
    cmp rsi, r11
    je .mr_check
.f_loop_body:
    cmp r12, [r14 + (psq_tab - base_primes) + rsi*8]
    jb .f_end               ; p^2 > remaining

    mov rax, r12
    imul rax, [r14 + (pinv_tab - base_primes) + rsi*8]
    cmp rax, [r14 + (plim_tab - base_primes) + rsi*8]
    ja .n_p                 ; p does not divide remaining

    ; p divides remaining; RAX is already remaining / p
    mov ebx, [r14 + rsi*4]  ; p
    mov r15, [r14 + (pinv_tab - base_primes) + rsi*8]
    mov r10, [r14 + (plim_tab - base_primes) + rsi*8]
    mov r8, 1
    mov r9, 1
.d_loop:
    mov r12, rax
    mov rax, r9
    mul rbx
    jc .ovf
    mov r9, rax             ; p^k
    add r8, r9              ; 1 + p + ... + p^k
    jc .ovf

    mov rax, r12
    imul rax, r15
    cmp rax, r10
    jbe .d_loop             ; still divisible: RAX = next quotient

    mov rax, r13
    mul r8
    jc .ovf
    mov r13, rax
    cmp r13, rbp
    ja .no_match            ; sigma already exceeds N + M: cannot match
    ; same divisibility + lower-bound abort as after the 2-strip
    mov rax, rbp
    xor edx, edx
    div r13
    test rdx, rdx
    jnz .no_match
    cmp r12, 1
    jbe .n_p
    lea rdx, [r12 + 1]
    cmp rax, rdx
    jb .no_match
    cmp r12, 0x100000       ; cofactor still > 2^20: re-arm the
    jbe .n_p                ; Miller-Rabin trigger 128 primes ahead
    lea r11, [rsi + 128]
.n_p:
    inc rsi
    jmp .f_loop
.f_end:
    cmp r12, 1
    jbe .done
    inc r12
    mov rax, r13
    mul r12
    jc .ovf
    mov r13, rax
.done:
    mov rax, r13
    sub rax, rcx

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    pop rbp
    ret
.no_match:
    xor eax, eax
    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    pop rbp
    ret
.ovf:
    mov rax, -1
    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    pop rbp
    ret

; Reached (once armed) after 128 consecutive base primes failed to divide
; the cofactor: it then very likely has no small factor at all, so settle
; it with one Miller-Rabin test instead of scanning on to sqrt(cofactor).
.mr_check:
    mov r11, -1             ; disarm (fires at most once per arming)
    cmp r12, 0x100000       ; small cofactor: rest of scan is cheap
    jbe .f_loop_body
    mov rax, r12
    call mr_is_prime        ; preserves everything except RAX
    test rax, rax
    jnz .f_end              ; prime: .f_end applies sigma *= (cofactor+1)
    jmp .f_loop_body        ; composite: resume scanning where we left off

; -------------------------------------------------------------------
; mr_is_prime - RAX = n (odd, 2^20 < n < 2^43) -> RAX = 1 if prime,
; else 0. Deterministic Miller-Rabin with bases {2,3,5,7,11,13,17},
; exact for all n < 341,550,071,728,321 (Jaeschke). The 5-base set used
; by the 2^32 build only reaches 2,152,302,898,747, which is below the
; largest cofactor here: sigma(N) can exceed 6N for superabundant N,
; so M = sigma(N) - N runs to roughly 5.7e12 at the top of the 2^40
; range. Arithmetic is Montgomery multiplication mod n (R = 2^64),
; reusing the Newton-iteration inverse trick from the reciprocal
; tables. Preserves all registers except RAX.
; -------------------------------------------------------------------
mr_is_prime:
    push rbx
    push rcx
    push rdx
    push rsi
    push rdi
    push r8
    push r9
    push r10
    push r11
    push r12
    push r13
    push r14
    push r15

    mov rbx, rax            ; n
    ; n-1 = d * 2^s (done first, while CL is still free for the shift)
    lea r8, [rbx - 1]
    bsf rcx, r8             ; s = tzcnt(n-1)
    mov r9, rcx
    shr r8, cl              ; d
    ; Newton: rcx = n^-1 mod 2^64 (5 doublings from seed n)
    mov rcx, rbx
    mov r10, 5
.newton:
    mov rax, rbx
    imul rax, rcx
    mov rdx, 2
    sub rdx, rax
    imul rcx, rdx           ; x *= 2 - n*x
    dec r10
    jnz .newton
    neg rcx                 ; -(n^-1) mod 2^64, as REDC wants
    ; rsi = R mod n (= mont(1)), rdi = n - rsi (= mont(n-1))
    xor edx, edx
    mov rax, -1
    div rbx
    lea rsi, [rdx + 1]      ; 2^64 mod n; n odd, so never wraps to 0
    mov rdi, rbx
    sub rdi, rsi
    ; bases 2,3,5,7,11,13,17 - one byte each, consumed low byte first
    mov r10, 0x110D0B07050302
.base_loop:
    movzx r11, r10b         ; a
    shr r10, 8
    mov rax, r11            ; mont(a) = a * (R mod n) mod n (fits: < 2^48)
    mul rsi
    div rbx                 ; RDX was 0 after the mul, so this is exact
    mov r11, rdx            ; base accumulator = mont(a)
    mov r12, rsi            ; x = mont(1)
    mov r13, r8             ; remaining bits of d
.pow_loop:
    test r13, 1
    jz .pow_sq
    mov rax, r12            ; x = REDC(x * base)
    mul r11
    mov r14, rax            ; lo
    mov r15, rdx            ; hi
    imul rax, rcx           ; m = lo * n'
    mul rbx                 ; RDX = high(m*n)
    neg r14                 ; CF = (lo != 0)
    adc r15, rdx            ; hi + high(m*n) + carry, in [0, 2n)
    mov r12, r15
    sub r15, rbx
    cmovnc r12, r15         ; subtract n if the sum was >= n
.pow_sq:
    shr r13, 1
    jz .pow_done
    mov rax, r11            ; base = REDC(base^2)
    mul r11
    mov r14, rax
    mov r15, rdx
    imul rax, rcx
    mul rbx
    neg r14
    adc r15, rdx
    mov r11, r15
    sub r15, rbx
    cmovnc r11, r15
    jmp .pow_loop
.pow_done:
    cmp r12, rsi            ; a^d == 1 ?
    je .base_pass
    cmp r12, rdi            ; a^d == n-1 ?
    je .base_pass
    mov r13, r9             ; else square s-1 more times,
    dec r13                 ; looking for n-1
.sq_loop:
    test r13, r13
    jz .composite
    mov rax, r12            ; x = REDC(x^2)
    mul r12
    mov r14, rax
    mov r15, rdx
    imul rax, rcx
    mul rbx
    neg r14
    adc r15, rdx
    mov r12, r15
    sub r15, rbx
    cmovnc r12, r15
    cmp r12, rdi
    je .base_pass
    dec r13
    jmp .sq_loop
.composite:
    xor eax, eax
    jmp .mr_ret
.base_pass:
    test r10, r10
    jnz .base_loop
    mov eax, 1
.mr_ret:
    pop r15
    pop r14
    pop r13
    pop r12
    pop r11
    pop r10
    pop r9
    pop r8
    pop rdi
    pop rsi
    pop rdx
    pop rcx
    pop rbx
    ret
