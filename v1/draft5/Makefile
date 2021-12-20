generic: generic.cpp
	g++ -std=c++20 -Wall generic.cpp -o generic

bug: bug.cpp
	g++ -std=c++20 -Wall -fconcepts-diagnostics-depth=3 bug.cpp -o bug 

bug2: bug2.cpp
	g++ -std=c++20 -Wall -fconcepts-diagnostics-depth=3 bug2.cpp -o bug2 

module: module.cpp test.cpp part1.cpp part2.cpp
	g++ -std=c++20 -Wall -fmodules-ts part1.cpp part2.cpp test.cpp module.cpp  -o module

# bug: bug.cpp
# 	clang++ -std=c++20 -Wall bug.cpp -o bug 
