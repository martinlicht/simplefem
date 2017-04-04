
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



indexrange.o: indexrange.hpp indexrange.cpp 
	$(CC) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CC) indexmap.cpp -c -o indexmap.o 

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CC) multiindex.cpp -c -o multiindex.o 

generateindexmaps.o: generateindexmaps.hpp generateindexmaps.cpp 
	$(CC) generateindexmaps.cpp -c -o generateindexmaps.o 

generatemultiindices.o: generatemultiindices.hpp generatemultiindices.cpp 
	$(CC) generatemultiindices.cpp -c -o generatemultiindices.o 

objects = indexrange.o indexmap.o multiindex.o generateindexmaps.o generatemultiindices.o

    
    

    
all: $(objects) 

include ../makefile.clean
