0000000000007cca <_ZN6Qubits6pauliXEh>:
    7cca:	53                   	push   %rbx
    7ccd:	49 89 f9             	mov    %rdi,%r9
    7cd0:	49 8b 41 10          	mov    0x10(%r9),%rax
    7cd4:	41 b8 10 00 00 00    	mov    $0x10,%r8d
    7cda:	89 f1                	mov    %esi,%ecx
    7cdc:	41 d3 e0             	shl    %cl,%r8d
    7cdf:	4d 63 c0             	movslq %r8d,%r8
    7ce2:	49 8b 51 10          	mov    0x10(%r9),%rdx
    7ce6:	4c 01 c2             	add    %r8,%rdx
    7ce9:	4d 89 c7             	mov    %r8,%r15
    7cec:	4d 01 c7             	add    %r8,%r15
    7cef:	4d 89 fa             	mov    %r15,%r10
    7cf2:	4d 03 51 10          	add    0x10(%r9),%r10
    7cf6:	4d 8b 59 08          	mov    0x8(%r9),%r11
    7cfa:	49 c1 e3 04          	shl    $0x4,%r11
    7cfe:	4d 03 59 10          	add    0x10(%r9),%r11
    7d02:	66 0f 28 12          	movapd (%rdx),%xmm2
    7d06:	66 0f 28 08          	movapd (%rax),%xmm1
    7d0a:	66 0f 29 10          	movapd %xmm2,(%rax)
    7d0e:	66 0f 29 0a          	movapd %xmm1,(%rdx)
    7d12:	48 83 c0 10          	add    $0x10,%rax
    7d16:	48 83 c2 10          	add    $0x10,%rdx
    7d1a:	4c 39 d2             	cmp    %r10,%rdx
    7d1d:	72 e3                	jb     7d02 <_ZN6Qubits6pauliXEh+0x38>
    7d1f:	4c 01 c2             	add    %r8,%rdx
    7d22:	4c 01 c0             	add    %r8,%rax
    7d25:	4d 01 fa             	add    %r15,%r10
    7d28:	4c 39 da             	cmp    %r11,%rdx
    7d2b:	72 d5                	jb     7d02 <_ZN6Qubits6pauliXEh+0x38>
    7d2f:	5b                   	pop    %rbx
    7d30:	c3                   	retq   
