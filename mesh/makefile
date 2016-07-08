
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

unittest.coordinates.out: coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp ../combinatorics/indexmap.o ../combinatorics/indexrange.o -o unittest.coordinates.out

simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp 
	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

generatesimplicialmesh.o: generatesimplicialmesh.hpp generatesimplicialmesh.cpp 
	$(CC) generatesimplicialmesh.cpp -c -o generatesimplicialmesh.o 

all: coordinates.o unittest.coordinates.out simplicialmesh.o generatesimplicialmesh.o

clean:
	rm *.exe *.o *.gch *.exe.stackdump *.out 