	// Implementation of Qubits::pauliX(std::uint8_t targ)
	//.global _ZN6Qubits6pauliXEh

	///\todo Find a systematic way of getting these mangled names
	//.global _ZN6QubitsIL4Type1EE6pauliXEh
	.global pauliX_x86

	#include "tools.S"
	#include "gate.S"

	//
	// Qubits::pauliX(std::uint8_t targ)
	//
	// The first parameter is 'this', which is a
	// pointer argument that ends up in rdi. The next
	// argument (targ) is stored in rsi.
	//
	// The implementation is as follows: rax and rdx
	// store the addresses of the elements which must
	// be swapped (which are k amplitudes apart,
	// meaning 16*k bytes. The factor of 16 is omitted
	// in the following logical addresses). This
	// distance is called the stride. The swaps are
	// arranged into a number of blocks as follows.
	// First, rax is incremented from 0 to k-1,
	// and rdx is incremented from k to 2*k-1. 
	//
	// When rdx reaches 2*k, the routine sets rax to
	// equal rdx. Then the new rdx is obtained by
	// adding k to the new rax (so that they are
	// separated by the stride k). This check
	// constitutes the outer loop. The inner loop
	// is obtained by incrementing rax and rdx by
	// 1 (until the new limit of 4*k is reached by
	// rdx).
	//
	// The loop exits when rax is set to an address
	// outside the state vector (i.e. it is equal
	// to the state vector size plus an offset)
	//
	///\todo Redo this comment to improve clarity
	//
//_ZN6Qubits6pauliXEh:
	//_ZN6QubitsIL4Type1EE6pauliXEh:
pauliX_x86:
//	print %rax
//	ret
	// Read state pointer and indices from Qubits class
	START_GATE

	// Main loop starts here
	START_LOOP

	// Swap the amplitudes and increment pointers
	movapd (%rax), %xmm1
	movapd (%rdx), %xmm2
	movapd %xmm2, (%rax)
	movapd %xmm1, (%rdx)
	
	// Loop overhead, clean-up and return
	END_GATE
