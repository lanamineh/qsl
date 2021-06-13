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
 * \file colours.hpp
 *
 * \brief Colours codes for printing to the termial
 *
 */

#ifndef COLOURS_HPP
#define COLOURS_HPP

/// Requires <string>, use of std::string

namespace Colour {

    const std::string GREEN{"\033[1;32m"};
    const std::string YELLOW{"\033[1;33m"};
    const std::string BLUE{"\033[1;34m"};
    const std::string CYAN{"\033[1;36m"};
    const std::string ORANGE{"\033[1;38;5;208m"};
    const std::string PINK{"\033[1;38;5;207m"};
    const std::string PALE{"\033[1;38;5;159m"};
    const std::string RESET{"\033[0m"};

    const std::string Address{CYAN};    
    const std::string Opcode{GREEN};
    const std::string ModRM{YELLOW};
    const std::string Register{YELLOW};
    const std::string RegInd{PINK};
    const std::string SIB{BLUE};
    const std::string Displacement{ORANGE};
    const std::string Immediate{PALE};
    const std::string Reset{RESET};
};

#endif
