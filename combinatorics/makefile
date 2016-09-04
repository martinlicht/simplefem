
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



indexrange.o: indexrange.hpp indexrange.cpp 
	$(CC) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CC) indexmap.cpp -c -o indexmap.o 

generateindexmaps.o: generateindexmaps.hpp generateindexmaps.cpp 
	$(CC) generateindexmaps.cpp -c -o generateindexmaps.o 

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CC) multiindex.cpp -c -o multiindex.o 

generatemultiindices.o: generatemultiindices.hpp generatemultiindices.cpp 
	$(CC) generatemultiindices.cpp -c -o generatemultiindices.o 

objects = indexrange.o indexmap.o generateindexmaps.o multiindex.o generatemultiindices.o

    
    
unittest.indexmap.out: indexmap.o indexrange.o indexmap.hpp unittest.indexmap.cpp 
	$(CC) unittest.indexmap.cpp indexmap.o indexrange.o -o unittest.indexmap.out

unittest.generateindexmaps.out: generateindexmaps.o indexmap.o indexrange.o indexmap.hpp unittest.generateindexmaps.cpp 
	$(CC) unittest.generateindexmaps.cpp generateindexmaps.o indexmap.o indexrange.o -o unittest.generateindexmaps.out

unittest.generatemultiindices.out: generatemultiindices.o indexmap.o indexrange.o multiindex.o unittest.generatemultiindices.cpp 
	$(CC) unittest.generatemultiindices.cpp generatemultiindices.o multiindex.o indexmap.o indexrange.o -o unittest.generatemultiindices.out

unittests = unittest.indexmap.out unittest.generateindexmaps.out unittest.generatemultiindices.out
    

    
all: $(objects) $(unittests) 

clean: 
	rm -f *.o *.exe *.exe.stackdump *.gch *.out