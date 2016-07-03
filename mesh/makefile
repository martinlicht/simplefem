
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp 
	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

unittest.coordinates.out: coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp ../combinatorics/indexmap.o ../combinatorics/indexrange.o -o unittest.coordinates.out

all: coordinates.o simplicialmesh.o unittest.coordinates.out
