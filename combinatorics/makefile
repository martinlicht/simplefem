
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

indexrange.o: indexrange.hpp indexrange.cpp 
	$(CC) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CC) indexmap.cpp -c -o indexmap.o 

unittest.indexmap.out: indexmap.o indexrange.o indexmap.hpp unittest.indexmap.cpp 
	$(CC) unittest.indexmap.cpp indexmap.o indexrange.o -o unittest.indexmap.out

generateindexmaps.o: generateindexmaps.hpp generateindexmaps.cpp 
	$(CC) generateindexmaps.cpp -c -o generateindexmaps.o 

unittest.generateindexmaps.out: generateindexmaps.o indexmap.o indexrange.o indexmap.hpp unittest.generateindexmaps.cpp 
	$(CC) unittest.generateindexmaps.cpp generateindexmaps.o indexmap.o indexrange.o -o unittest.generateindexmaps.out

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CC) multiindex.cpp -c -o multiindex.o 

all: indexrange.o indexmap.o unittest.indexmap.out multiindex.o generateindexmaps.o unittest.generateindexmaps.out

clean: 
	rm *.o *.exe *.exe.stackdump *.gch *.out