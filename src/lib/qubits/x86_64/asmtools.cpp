/* 
 * Copyright (C) 2020 Lana Mineh and John Scott.
 *
 * This file is part of QSL, the quantum computer simulator.
 *
 * QSL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * QSL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with QSL.  If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * \file asmtools.cpp
 * \brief Tools for helping program assembly functions
 *
 */

#include <iostream>

extern "C" void finish() {
    while(true);
}

extern "C" void print(std::size_t val) {
    std::cout << std::hex << val << std::endl;
}

extern "C" void print_double(double val) {
    std::cout << std::dec << val << std::endl;
}
