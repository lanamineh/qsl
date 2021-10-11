	// Implementation of Qubits::pauliX(std::uint8_t targ)
	//.global _ZN6Qubits7rotateXEhd

	///\todo Work out mangled name
	//.global _ZN6QubitsIL4Type1EE7rotateXEhd
	.global rotateX_x86
	
	#include "gate.S"
	#include "tools.S"
	
	// This macro reads the rotation angle from the
	// argument list
	.macro GET_ANGLE

	// Store the integer constant 0.5 in xmm1 (see Agner Fog manual)
	pcmpeqw %xmm1,%xmm1
	psllq $55,%xmm1
	psrlq $2,%xmm1
	
	// Get the angle argument (from %xmm0), divide by 2
	mulsd %xmm0, %xmm1

	// In the following, it is important to put the
	// l suffix on the fld and fstp instructions so
	// that the full 64 bit value is copied to and
	// from the x87 stack.
	sub $24, %rsp
	movsd %xmm1, (%rsp) 	// Move xmm1 to the stack
	fldl (%rsp)   		// st0 = xmm1
	fld %st(0)		// Duplicate stack top

	fcos			// Replace st0 with cos(st0)
	fstpl (%rsp)     	// Store cos back on stack, pop x87 stack
	fwait
	movsd (%rsp), %xmm4 	// Move cos value to xmm5	
	
	fsin			// Replace st1 (now st0) with sin(st1)
	fstpl (%rsp)     	// Store sin back on stack, pop the stack
	movsd (%rsp), %xmm5 	// Move sin value to xmm6	
	
	add $24, %rsp
	
	// End macro
	.endm
	
	//
	// Qubits::rotateX(std::uint8_t targ, double angle)
	//
	// This routine uses the same method as pauliX to
	// access the state vector entries. It takes
	// slightly longer to run due to the arithmetic
	// inside the inner loop. There is more optimisation
	// to be done (e.g. get rid of vex instructions) but
	// it is already a bit faster than the gcc version.
	//
//_ZN6Qubits7rotateXEhd:
	//_ZN6QubitsIL4Type1EE7rotateXEhd:
rotateX_x86:

	// Read the state vector pointer and indices from Qubits class
	START_GATE

	// Get the rotation angle
	GET_ANGLE
	
	// Put the content of the gate after this line
	START_LOOP

	// 4,6 = cos, 5,7 = sin
	movsd (%rax),%xmm2    // Real part a0 -> 2
	movsd (%rdx),%xmm1    // Real part a1 -> 1
	movsd 0x8(%rdx),%xmm0 // Imag part a1 -> 0
	movsd 0x8(%rax),%xmm3 // Imag part a0 -> 3

	movapd %xmm4,%xmm6 // Copy cos to 6
	mulsd %xmm2,%xmm6 // 6 = a0.real * cos
	mulsd %xmm5,%xmm2 // 2 = a0.real * sin
	movapd %xmm0,%xmm7 // 7 = a1.imag
	mulsd %xmm5,%xmm7 // 7 = a1.imag * sin
	mulsd %xmm4,%xmm0 // 0 = a1.imag * cos
	addsd %xmm7,%xmm6 // 6 = a0.real * cos + a1.imag * sin (new a0.real) 

	movapd %xmm5,%xmm7 // Copy sin to 7
	mulsd %xmm1,%xmm7 // 7 = a1.real * sin
	subsd %xmm2,%xmm0 // 0 = a1.imag * cos - a0.real * sin (new a1.imag)
	mulsd %xmm4,%xmm1 // 1 = a1.real * cos

	movsd %xmm6,(%rax) // Move 6 to a0.real
	
	movapd %xmm3,%xmm6 // 6 = a0.imag
	mulsd %xmm4,%xmm6 // 6 = a0.imag * cos
	mulsd %xmm5,%xmm3 // 3 = a0.imag * sin
	subsd %xmm7,%xmm6 // 6 = a0.imag * cos - a1.real * sin (new a0.imag)

	addsd %xmm3,%xmm1 // 1 = a0.imag * sin + a1.real * cos (a1.real)
	
	movsd %xmm6, 0x8(%rax) // Store a0.imag
	movsd %xmm1, (%rdx) // Store a1.real
	movsd %xmm0, 0x8(%rdx) // Store a1.imag
	
	// Loop overhead, clean up, return
	END_GATE
