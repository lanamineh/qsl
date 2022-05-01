CC=g++-10

test.o: check
	${CC} -c -std=c++20 -Wall -Wextra test.cpp -o test.o

check:
	cppcheck --enable=all --suppress=missingIncludeSystem qsl.hpp test.cpp

docs:
	doxygen
