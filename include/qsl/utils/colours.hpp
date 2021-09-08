/*
 *  Authors: Lana Mineh and John Scott
 *  Copyright 2021 Phasecraft Ltd. and John Scott
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
 * \file colours.hpp
 *
 * \brief Colours codes for printing to the termial
 *
 */

#ifndef QSL_COLOURS_HPP
#define QSL_COLOURS_HPP

/// Requires <string>, use of std::string
namespace qsl {
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
}
    
#endif
