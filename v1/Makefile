all:
	g++-10 -std=c++20 qsl.hpp -Wall -Wextra

check:
	cppcheck --enable=all --suppress=missingIncludeSystem qsl.hpp
