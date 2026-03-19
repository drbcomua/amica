; -------------------------------------------------------------------
; amica_fast_linux.asm - Multi-threaded Segmented Sieve for Linux x64
; Supports full 32-bit unsigned range (0 to 4,294,967,295)
; Optimized for 16-thread Intel Core i5-14400F (AVX2)
; Integrates Odds-Only Prime Sieve & Segmented Divisor Sieve
; Added: Real-time Prime Factorization & Erdős Type Classification
; -------------------------------------------------------------------
default rel

extern printf
extern sprintf
extern fflush
extern pthread_create
extern pthread_join
extern exit

section .data
    limit          equ 4294967295
    sqrt_limit     equ 262144      ; Cover factors up to 68.7 billion for S(N)
    segment_size   equ 131072      ; 128 KB (Fits in L2 cache)
    segment_bits   equ segment_size * 8
    num_threads    equ 16          ; Matches Core i5-14400F logical cores
    
    ; Output Formatters
    fmt            db "Primes found: %llu, Largest: %u", 10, 0
    overflow_fmt   db "Lower int to start in 128-bit range: %llu", 10, 0
    summary_fmt    db "Checked pairs in the range 2..4294967295; found: %lld amicable pairs", 10, 0
    
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
    
    print_lock     dd 0

section .bss
    alignb 32
    base_primes       resd 30000      
    base_prime_cnt    resd 1
    
    ; Thread management
    thread_handles    resq num_threads
    
    ; Statically allocated memory
    alignb 32
    segment_buffers   resb num_threads * segment_size
    alignb 32
    next_indices      resq num_threads * 30000
    alignb 32
    base_sieve_buffer resb 32768      
    
    ; Segmented Divisor Sieve Buffers (8MB rem + 16MB sigma per thread)
    alignb 32
    chunk_rem         resd 33554432   
    alignb 32
    chunk_sigma       resq 33554432   
    
    ; Atomic Dynamic Work Queue & Synchronization
    alignb 8
    global_task_start resq 1
    global_prime_cnt  resq 1
    global_pair_cnt   resq 1
    global_largest    resd 1
    overflow_reported resq 1

section .text
    global main

main:
    push rbp
    mov rbp, rsp
    sub rsp, 80             ; Keep stack 16-byte aligned

    mov dword [rel print_lock], 0

    ; 1. Setup base sieve (0 to 262,144)
    lea rdi, [rel base_sieve_buffer]
    mov rcx, 4096           
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
    cmp rsi, 512            
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

    ; Initialize Globals
    mov qword [rel global_task_start], 0
    mov qword [rel global_prime_cnt], 0
    mov dword [rel global_largest], 0
    mov qword [rel global_pair_cnt], 0
    mov qword [rel overflow_reported], 0

    ; -------------------------------------------------------------------
    ; WAVE 1: Prime Sieve
    ; -------------------------------------------------------------------
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
    mov edx, [rel global_largest]
    xor eax, eax
    call printf

    ; -------------------------------------------------------------------
    ; WAVE 2: Amicable Pair Search (Segmented Divisor Sieve)
    ; -------------------------------------------------------------------
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
    shl rax, 17             
    lea r14, [rel segment_buffers]
    add r14, rax

    mov rax, rbp
    imul rax, 240000        
    lea r15, [rel next_indices]
    add r15, rax

    xor r12, r12            
    xor r13, r13            

.get_task:
    mov rax, 2097152
    lea rcx, [rel global_task_start]
    lock xadd [rcx], rax
    mov rsi, rax            
    
    mov rax, limit
    cmp rsi, rax
    jae .thread_end         

    mov r10, rsi
    add r10, 2097152
    cmp r10, rax
    jbe .init_indices
    mov r10, rax
    inc r10                 

.init_indices:
    xor rbx, rbx
    mov ecx, dword [rel base_prime_cnt]
    lea rdx, [rel base_primes]
.idx_loop:
    mov eax, [rdx + rbx*4]  
    mov r8, rax
    mul rax                 
    cmp rax, rsi
    jae .store_idx
    mov rax, rsi
    xor rdx, rdx
    div r8
    mul r8                  
    cmp rax, rsi
    jae .chk_odd_idx
    add rax, r8
.chk_odd_idx:
    test rax, 1
    jnz .store_idx
    add rax, r8
.store_idx:
    mov [r15 + rbx*8], rax
    inc rbx
    lea rdx, [rel base_primes]
    cmp ebx, ecx
    jb .idx_loop

    mov rdi, r14
    mov rcx, 16384          
    mov rax, 0xFFFFFFFFFFFFFFFF
    rep stosq

    test rsi, rsi
    jnz .sieve
    btr qword [r14], 0      
    inc r12                 

.sieve:
    xor rbx, rbx
    mov ecx, dword [rel base_prime_cnt]
    lea r8, [rel clear_masks]
    lea r11, [rel base_primes]
.prime_loop:
    mov r9, [r15 + rbx*8]
    mov eax, [r11 + rbx*4]
    
    cmp r9, r10
    jae .skip_p
    lea rax, [rax + rax]       
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
    mov rcx, 4096
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

    ; Find Largest Prime in Chunk
    lea rdi, [r14 + 131071]
    mov rcx, 131072
.f_last:
    movzx eax, byte [rdi]
    test eax, eax
    jnz .found_p
    dec rdi
    dec rcx
    jnz .f_last
    jmp .get_task
.found_p:
    bsr eax, eax
    dec rcx
    shl rcx, 4              
    lea rax, [rax + rax]
    add rax, rcx
    add rax, rsi
    inc rax                 
    cmp rax, r10
    jae .retry_last
    cmp eax, r13d
    cmova r13d, eax
    jmp .get_task
.retry_last:
    dec rdi
    dec rcx
    jnz .f_last

.thread_end:
    lea rcx, [rel global_prime_cnt]
    mov rax, r12
    lock xadd [rcx], rax
    
.upd_largest:
    mov eax, dword [rel global_largest]
    cmp r13d, eax
    jbe .done_end
    mov ecx, eax
    mov edx, r13d
    lock cmpxchg dword [rel global_largest], edx
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
    shl rax, 21             
    
    lea r14, [rel chunk_rem]
    lea r8, [rax * 4]
    add r14, r8             
    
    lea r15, [rel chunk_sigma]
    lea r8, [rax * 8]
    add r15, r8             

    xor r12, r12            

.get_task:
    mov rax, 2097152
    lea rcx, [rel global_task_start]
    lock xadd [rcx], rax
    mov rsi, rax            
    
    mov rax, limit
    cmp rsi, rax
    jae .thread_end         

    mov r10, rsi
    add r10, 2097152
    cmp r10, rax
    jbe .set_end
    mov r10, rax
    inc r10                 
.set_end:
    mov rax, r10
    sub rax, rsi
    mov [rsp + 64], rax     

    ; 1. Initialize Sigma Array to 1
    mov rdi, r15
    mov rcx, [rsp + 64]
    mov rax, 1
    rep stosq

    ; 2. Initialize Remainder Array to N
    mov rdi, r14
    mov rcx, [rsp + 64]
    mov r8d, esi            
.init_rem:
    mov [rdi], r8d
    add rdi, 4
    inc r8d
    dec rcx
    jnz .init_rem

    ; 3. Segmented Divisor Sieve
    xor r8d, r8d            
.sieve_2_loop:
    cmp r8d, dword [rsp + 64]
    jae .sieve_odd_primes
    
    mov eax, [r14 + r8*4]   
    test eax, eax
    jz .next_2              
    
    mov r10, 2              
    mov r11, 3              
    shr eax, 1              
    mov ecx, eax
    
    test eax, 1
    jnz .end_div_2
    
.div_2_while:
    shr eax, 1
    mov ecx, eax
    shl r10, 1
    add r11, r10
    test eax, 1
    jz .div_2_while
    
.end_div_2:
    mov [r14 + r8*4], ecx
    mov rax, [r15 + r8*8]
    mul r11
    mov [r15 + r8*8], rax
    
.next_2:
    add r8d, 2
    jmp .sieve_2_loop

.sieve_odd_primes:
    mov ebx, dword [rel base_prime_cnt]
    xor r13d, r13d          
    lea rdi, [rel base_primes]
.prime_loop:
    cmp r13d, ebx
    jae .primes_done
    
    mov r9d, [rdi + r13*4]  
    cmp r9d, 65536
    jae .primes_done        

    mov eax, esi
    xor edx, edx
    div r9d
    test edx, edx
    jz .mod_zero
    mov eax, r9d
    sub eax, edx
    jmp .set_start
.mod_zero:
    xor eax, eax
.set_start:
    mov r8d, eax            
    
    test esi, esi
    jnz .sieve_p_loop
    mov r8d, r9d            
    
.sieve_p_loop:
    cmp r8d, dword [rsp + 64]
    jae .next_prime
    
    mov eax, [r14 + r8*4]   
    test eax, eax
    jz .next_i
    
    xor edx, edx
    div r9d                 
    
    mov r10, r9             
    lea r11, [r9 + 1]       
    mov ecx, eax            
    
    test eax, eax
    jz .end_div_while
    
.div_while:
    xor edx, edx
    div r9d
    test edx, edx
    jnz .end_div_while      
    
    mov ecx, eax
    imul r10, r9            
    add r11, r10            
    jmp .div_while
    
.end_div_while:
    mov [r14 + r8*4], ecx   
    mov rax, [r15 + r8*8]
    mul r11
    mov [r15 + r8*8], rax   
    
.next_i:
    add r8d, r9d
    jmp .sieve_p_loop

.next_prime:
    inc r13d
    jmp .prime_loop

.primes_done:

    ; 4. Calculate S(N)
    xor r8d, r8d            
.finish_loop:
    cmp r8d, dword [rsp + 64]
    jae .am_scan
    
    mov eax, esi
    add eax, r8d            
    cmp eax, 1
    jbe .skip_finish
    
    mov ecx, [r14 + r8*4]   
    mov r11, [r15 + r8*8]   
    
    cmp ecx, 1
    jbe .no_rem
    mov eax, ecx
    inc rax
    mul r11                 
    mov r11, rax
.no_rem:
    
    mov eax, esi
    add eax, r8d
    sub r11, rax            
    mov [r15 + r8*8], r11   
    
.skip_finish:
    inc r8d
    jmp .finish_loop

    ; 5. Instant Memory Lookup Amicable Scan
.am_scan:
    xor rdi, rdi            
.am_scan_loop:
    cmp rdi, [rsp + 64]
    jae .get_task
    
    mov eax, esi
    add eax, edi            
    cmp eax, 1
    jbe .skip_am
    
    mov r13, [r15 + rdi*8]  
    cmp r13, rax            
    jbe .skip_am
    
    mov r9, r13
    sub r9, rsi
    cmp r9, [rsp + 64]
    jae .compute_ssn
    
    mov r10, [r15 + r9*8]   
    jmp .check_pair
    
.compute_ssn:
    mov rcx, r13            
    call SumProperDivisors
    cmp rax, -1
    je .handle_ovf
    mov r10, rax
    
.check_pair:
    mov eax, esi
    add eax, edi
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
    
    mov eax, esi
    add eax, edi
    push rdi
    push rsi
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
; PrintAmicablePair - Factorizes and formats pairs with Erdős (i,j)
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
    ; Format into [rbp-1320]
    lea rdi, [rbp-1320]     
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
    
    lea r12, [rbp-1320]
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
    lea rsi, [rbp-1320]
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
    mov r11, 3      ; p = 3
.loop3:
    mov rax, r11
    mul r11
    cmp rax, r8
    ja .done3
    
    mov rax, r8
    xor rdx, rdx
    div r11
    test rdx, rdx
    jnz .next_p
    
    xor r10, r10    ; exp = 0
.div_p:
    mov rax, r8
    xor rdx, rdx
    div r11
    test rdx, rdx
    jnz .done_p
    inc r10
    mov r8, rax     ; N = N / p
    jmp .div_p
.done_p:
    mov qword [rsi + r9*8], r11
    mov dword [r12 + r9*4], r10d
    inc r9
    
.next_p:
    add r11, 2
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
; SumProperDivisors - Fallback Trial Division
; -------------------------------------------------------------------
SumProperDivisors:
    push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15              

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
    
.skip_2:
    xor rsi, rsi
    lea r14, [rel base_primes]  
    mov edi, dword [rel base_prime_cnt]
    
.f_loop:
    cmp esi, edi
    jae .ovf
    mov ebx, [r14 + rsi*4]
    
    mov rax, rbx
    mul rbx
    cmp rax, r12
    ja .f_end
    
    mov rax, r12
    shr rax, 32
    test rax, rax
    jnz .div64
    
    mov eax, r12d
    xor edx, edx
    div ebx
    jmp .div_check
    
.div64:
    mov rax, r12
    xor rdx, rdx
    div rbx
    
.div_check:
    test rdx, rdx
    jnz .n_p
    
    mov r8, 1
    mov r9, 1
.d_loop:
    mov r12, rax
    mov rax, r9
    mul rbx
    jc .ovf
    mov r9, rax
    add r8, r9
    jc .ovf
    
    mov rax, r12
    shr rax, 32
    test rax, rax
    jnz .d_div64
    
    mov eax, r12d
    xor edx, edx
    div ebx
    jmp .d_div_check
    
.d_div64:
    mov rax, r12
    xor rdx, rdx
    div rbx
    
.d_div_check:
    test rdx, rdx
    jz .d_loop
    
    mov rax, r13
    mul r8
    jc .ovf
    mov r13, rax
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
    ret