test.o: check
	g++-10 -c -std=c++20 -Wall -Wextra test.cpp -o test.o

check:
	cppcheck --enable=all --suppress=missingIncludeSystem qsl.hpp test.cpp
