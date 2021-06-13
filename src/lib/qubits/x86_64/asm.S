0000000000007bca <_ZN6Qubits6pauliXEh>:
7bca:	48 8b 47 10          	mov    0x10(%rdi),%rax
7bce:	49 c7 c0 10 00 00 00 	mov    $0x10,%r8
7bd5:	48 89 f1             	mov    %rsi,%rcx
7bd8:	49 d3 e0             	shl    %cl,%r8
7bdb:	48 8b 57 10          	mov    0x10(%rdi),%rdx
7bdf:	4c 01 c2             	add    %r8,%rdx
7be2:	4f 8d 14 00          	lea    (%r8,%r8,1),%r10
7be6:	4c 03 57 10          	add    0x10(%rdi),%r10
7bea:	4c 8b 5f 08          	mov    0x8(%rdi),%r11
7bee:	49 c1 e3 04          	shl    $0x4,%r11
7bf2:	4c 03 5f 10          	add    0x10(%rdi),%r11
	
	// Inner and outer loop start
7bf6:	66 0f 28 08          	movapd (%rax),%xmm1
7bfa:	66 0f 28 12          	movapd (%rdx),%xmm2
7bfe:	66 0f 29 10          	movapd %xmm2,(%rax)
7c02:	66 0f 29 0a          	movapd %xmm1,(%rdx)
7c06:	48 83 c0 10          	add    $0x10,%rax
7c0a:	48 83 c2 10          	add    $0x10,%rdx
7c0e:	4c 39 d2             	cmp    %r10,%rdx
7c11:	72 e3                	jb     7bf6 <_ZN6Qubits6pauliXEh+0x2c>
	// Inner loop end

7c13:	4c 01 c2             	add    %r8,%rdx
7c16:	4f 8d 14 42          	lea    (%r10,%r8,2),%r10
7c1a:	4c 01 c0             	add    %r8,%rax
7c1d:	4c 39 da             	cmp    %r11,%rdx
7c20:	72 d4                	jb     7bf6 <_ZN6Qubits6pauliXEh+0x2c>
	// outer loop end

7c22:	c3                   	retq   
