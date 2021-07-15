/**
 * \file test-utils.cpp
 * \brief Implementation of test utilities
 *
 */

/// Gets the bit in the nth position of val
unsigned getBit(unsigned val, unsigned n)
{
    return 1 & (val >> n);
}

/// Set the nth bit of val to b
unsigned setBit(unsigned val, unsigned n, unsigned b)
{
    if (b == 0) {
	val &= ~(1 << n);
    } else {
	val |= (1 << n);
    }
    return val;
}
