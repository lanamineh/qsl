LFLAGS = -L../build/lib/ -Wl,-rpath=../build/lib/
IFLAGS = -I../include/

basics: basics.cpp
	g++ -std=c++17 ${IFLAGS} ${LFLAGS} basics.cpp -o basics -lqsl
