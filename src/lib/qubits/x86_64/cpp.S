00000000000079e0 <_ZN6Qubits6pauliXEh>:
79e0:	4c 8b 57 08          	mov    0x8(%rdi),%r10
79e4:	41 b8 01 00 00 00    	mov    $0x1,%r8d
79ea:	89 f1                	mov    %esi,%ecx
79ec:	41 d3 e0             	shl    %cl,%r8d
79ef:	4d 85 d2             	test   %r10,%r10
79f2:	74 71                	je     7a65 <_ZN6Qubits6pauliXEh+0x85>
79f4:	4d 63 c0             	movslq %r8d,%r8
79f7:	49 89 f9             	mov    %rdi,%r9
79fa:	53                   	push   %rbx
79fb:	31 f6                	xor    %esi,%esi
79fd:	4d 89 c3             	mov    %r8,%r11
7a00:	4c 89 c7             	mov    %r8,%rdi
7a03:	4b 8d 1c 00          	lea    (%r8,%r8,1),%rbx
7a07:	49 c1 e3 05          	shl    $0x5,%r11
7a0b:	48 c1 e7 04          	shl    $0x4,%rdi
7a0f:	90                   	nop

	// Start of outer loop
7a10:	4d 85 c0             	test   %r8,%r8
7a13:	74 43                	je     7a58 <_ZN6Qubits6pauliXEh+0x78>
7a15:	49 8b 51 10          	mov    0x10(%r9),%rdx
7a19:	48 89 f0             	mov    %rsi,%rax
7a1c:	48 c1 e0 04          	shl    $0x4,%rax
7a20:	49 03 41 10          	add    0x10(%r9),%rax
7a24:	48 01 fa             	add    %rdi,%rdx
7a27:	48 89 d1             	mov    %rdx,%rcx
7a2a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)

	// Start of inner loop
7a30:	f3 0f 6f 12          	movdqu (%rdx),%xmm2
7a34:	f2 0f 10 08          	movsd  (%rax),%xmm1
7a38:	48 83 c0 10          	add    $0x10,%rax
7a3c:	48 83 c2 10          	add    $0x10,%rdx
7a40:	f2 0f 10 40 f8       	movsd  -0x8(%rax),%xmm0
7a45:	0f 11 50 f0          	movups %xmm2,-0x10(%rax)
7a49:	f2 0f 11 4a f0       	movsd  %xmm1,-0x10(%rdx)
7a4e:	f2 0f 11 42 f8       	movsd  %xmm0,-0x8(%rdx)
7a53:	48 39 c8             	cmp    %rcx,%rax
7a56:	75 d8                	jne    7a30 <_ZN6Qubits6pauliXEh+0x50>
	// End of inner loop
	
7a58:	48 01 de             	add    %rbx,%rsi
7a5b:	4c 01 df             	add    %r11,%rdi
7a5e:	4c 39 d6             	cmp    %r10,%rsi
7a61:	72 ad                	jb     7a10 <_ZN6Qubits6pauliXEh+0x30>
	// End of outer loop
	
7a63:	5b                   	pop    %rbx
7a64:	c3                   	retq   
7a65:	c3                   	retq   
7a66:	66 2e 0f 1f 84 00 00 	nopw   %cs:0x0(%rax,%rax,1)
7a6d:	00 00 00 
	
