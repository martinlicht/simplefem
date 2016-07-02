
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

indexrange.o: indexrange.hpp indexrange.cpp 
	$(CC) indexrange.cpp -c -o indexrange.o 

natmap.o: natmap.hpp natmap.cpp 
	$(CC) natmap.cpp indexrange.o -c -o natmap.o 

unittest_natmap: natmap.o indexrange.o natmap.hpp unittest_natmap.cpp 
	$(CC) unittest_natmap.cpp natmap.o indexrange.o -o unittest_natmap

all: indexrange.o natmap.o unittest_natmap
