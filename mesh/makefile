
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

unittest.coordinates.out: coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp coordinates.o ../combinatorics/*.o ../operators/*.o -o unittest.coordinates.out

simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp coordinates.o
	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

generatesimplicialmesh.o: generatesimplicialmesh.hpp generatesimplicialmesh.cpp simplicialmesh.o coordinates.o
	$(CC) generatesimplicialmesh.cpp -c -o generatesimplicialmesh.o 

unittest.simplicialmesh.out: simplicialmesh.o simplicialmesh.hpp unittest.simplicialmesh.cpp generatesimplicialmesh.o
	$(CC) unittest.simplicialmesh.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o ../operators/*.o ../combinatorics/*.o -o unittest.simplicialmesh.out

all: coordinates.o unittest.coordinates.out simplicialmesh.o unittest.simplicialmesh.out generatesimplicialmesh.o

clean:
	rm *.exe *.o *.gch *.exe.stackdump *.out 