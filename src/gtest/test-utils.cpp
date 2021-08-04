/*
 *  Copyright 2021 Lana Mineh and John Scott
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */
 
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
