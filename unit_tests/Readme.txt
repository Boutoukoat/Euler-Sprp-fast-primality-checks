

Unit tests : the corner-stone of good principles in software

./m_reg_test : check the assembly language trickeries, where 1 instruction of the x86-64 processor replaces a few lines a C code, and uses the fastest implementation. 

	Note that the preferred implementation is based on slowest cycle counts from https://www.agner.org/optimize/instruction_tables.pdf

	By toggling INLINE_ASM from 1 to 0, in m_reg.h, the corresponding C code is tested. This C code is valid on ARM, and hopefully on GPUs.

./generic_test : check the modular square algorithms as involved in sprp and fermat tests. 

	Generic code is a good compromize between the optimized versions and the proven-correct slow version.

