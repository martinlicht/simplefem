
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

natmap.o: natmap.hpp natmap.cpp 
	$(CC) natmap.cpp -c -o natmap.o 

unittest_natmap: natmap.o natmap.hpp unittest_natmap.cpp 
	$(CC) unittest_natmap.cpp natmap.o -o unittest_natmap

all: natmap.o unittest_natmap
