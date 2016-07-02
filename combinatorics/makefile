
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

indexrange.o: indexrange.hpp indexrange.cpp 
	$(CC) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CC) indexmap.cpp indexrange.o -c -o indexmap.o 

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CC) multiindex.cpp -c -o multiindex.o 

unittest.indexmap.out: indexmap.o indexrange.o indexmap.hpp unittest.indexmap.cpp 
	$(CC) unittest.indexmap.cpp indexmap.o indexrange.o -o unittest.indexmap.out

all: indexrange.o indexmap.o unittest.indexmap.out multiindex.o

clean: *.o *.exe *.exe.stackdump *.gch *.out 
	rm *.o *.exe *.exe.stackdump *.gch *.out