# Disassembly ASM versions
gdb -batch -ex "disassemble _ZN6QubitsIL4Type1EE6pauliXEh" "main.bin" > disasm/pauliX_asm.S
gdb -batch -ex "disassemble _ZN6QubitsIL4Type1EE7rotateXEhd" "main.bin" > disasm/rotateX_asm.S

# Disassemble compiler versions
gdb -batch -ex "disassemble _ZN6QubitsIL4Type0EE6pauliXEh" "main.bin" > disasm/pauliX_g++.S
gdb -batch -ex "disassemble _ZN6QubitsIL4Type0EE7rotateXEhd" "main.bin" > disasm/rotateX_g++.S

