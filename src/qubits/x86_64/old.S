0000000000007cca <_ZN6Qubits6pauliXEh>:
    7cca:	53                   	push   %rbx
    7ccb:	49 89 f9             	mov    %rdi,%r9
    7cce:	49 8b 41 10          	mov    0x10(%r9),%rax
    7cd2:	41 b8 10 00 00 00    	mov    $0x10,%r8d
    7cd8:	89 f1                	mov    %esi,%ecx
    7cda:	41 d3 e0             	shl    %cl,%r8d
    7cdd:	4d 63 c0             	movslq %r8d,%r8
    7ce0:	49 8b 51 10          	mov    0x10(%r9),%rdx
    7ce4:	4c 01 c2             	add    %r8,%rdx
    7ce7:	4d 89 c7             	mov    %r8,%r15
    7cea:	4d 01 c7             	add    %r8,%r15
    7ced:	4d 89 fa             	mov    %r15,%r10
    7cf0:	4d 03 51 10          	add    0x10(%r9),%r10
    7cf4:	4d 8b 59 08          	mov    0x8(%r9),%r11
    7cf8:	49 c1 e3 04          	shl    $0x4,%r11
    7cfc:	4d 03 59 10          	add    0x10(%r9),%r11
    7d00:	66 0f 28 12          	movapd (%rdx),%xmm2
    7d04:	66 0f 28 08          	movapd (%rax),%xmm1
    7d08:	66 0f 29 10          	movapd %xmm2,(%rax)
    7d0c:	66 0f 29 0a          	movapd %xmm1,(%rdx)
    7d10:	48 83 c0 10          	add    $0x10,%rax
    7d14:	48 83 c2 10          	add    $0x10,%rdx
    7d18:	4c 39 d2             	cmp    %r10,%rdx
    7d1b:	72 e3                	jb     7d00 <_ZN6Qubits6pauliXEh+0x36>
    7d1d:	4c 01 c2             	add    %r8,%rdx
    7d20:	4c 01 c0             	add    %r8,%rax
    7d23:	4d 01 fa             	add    %r15,%r10
    7d26:	4c 39 da             	cmp    %r11,%rdx
    7d29:	72 d5                	jb     7d00 <_ZN6Qubits6pauliXEh+0x36>
    7d2b:	5b                   	pop    %rbx
    7d2c:	c3                   	retq   
