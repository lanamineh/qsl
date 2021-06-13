	// Macro to print reg
	.macro print reg=%rax	// thing

	// Preserve registers
	push %rax
	push %rcx
	push %rdx
	push %rsi
	push %rdi
	push %r8
	push %r9
	push %r10
	push %r11

	// Call the function
	mov \reg, %rdi
	call print

	// Restore registers
	pop %r11
	pop %r10
	pop %r9
	pop %r8
	pop %rdi
	pop %rsi
	pop %rdx
	pop %rcx
	pop %rax
	
        .endm

	// Macro to print reg
	.macro print_double reg=%xmm0

	// Preserve registers
	push %rax
	push %rcx
	push %rdx
	push %rsi
	push %rdi
	push %r8
	push %r9
	push %r10
	push %r11
	
	// Call print function
	movsd \reg, %xmm0
	call print_double

	// Restore registers
	pop %r11
	pop %r10
	pop %r9
	pop %r8
	pop %rdi
	pop %rsi
	pop %rdx
	pop %rcx
	pop %rax
	
	
	.endm
	
