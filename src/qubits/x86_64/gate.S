	//
	// \file gate.S
	// \brief Contains the generic implementation for single qubit gates
	//
	//
	
	// Index setup for single qubit gate
	//
	// This macro reads the size of the state vector,
	// the state vector pointer, and various other
	// necessary index date from the Qubits object,
	// and stores them in the following registers:
	//
	//
	//
	.macro START_GATE

	// Get the pointer to the state vector,
	// the value of k (which gives the stride
	// length), the starting value for the
	// current limit, and the state vector size
	// (which is the upper limit)

	mov 0x10(%rdi), %rax  	// Store pointer to state vector

	// Compute 16*k = 16 << targ
	mov $0x10,%r8 		// Store literal 16
	mov %rsi,%rcx 		// Copy targ to ecx (for cl below)
	shl %cl,%r8 		// Shift 16 left by targ (stored in cl)
	
	// Work out initial value of rdx (state pointer + 16*k)
	mov 0x10(%rdi),%rdx
	add %r8, %rdx
	
	// Work out initial limit for rdx (state_ptr + 32*k), store in r10
	lea (%r8,%r8), %r10 // compute 2*r8, store in r10 
	add 0x10(%rdi), %r10
	
	// Work out limit for rax (state_ptr + state_size), store in r11
	mov 0x8(%rdi),%r11
	shl $0x4, %r11  
	add 0x10(%rdi),%r11

	// End macro
	.endm

	// This label is for the start of the loop
	.macro START_LOOP
0:
	// End macro
	.endm

	// This macro handles the jump
	// back to the top to repeat the inner loop
	//
	.macro END_GATE

	// Increment pointers
	add $0x10,%rax
	add $0x10,%rdx
	
	cmp %r10, %rdx 	// Check if rdx is still below current block limit
	jb 0b		// If it is, go back to 0
	// End of inner loop

	add %r8, %rdx 	// Set rdx = rax + 16*k (increment old rdx by 16*k)
	lea (%r10, %r8, 2), %r10 // Update the current limit in r10 by adding 32*k
	add %r8, %rax	// Increment rax to next block (by adding 16*k)
	
	cmp %r11, %rdx	// Check if rdx is below state vector limit
	jb 0b		// If it is, go back to 0
	// End of outer loop

	// Clean up and return
	ret

	// End macro
	.endm
