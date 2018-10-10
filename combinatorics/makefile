
# CXX = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



indexrange.o: indexrange.hpp indexrange.cpp 
	$(CXX) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CXX) indexmap.cpp -c -o indexmap.o 

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CXX) multiindex.cpp -c -o multiindex.o 

generateindexmaps.o: generateindexmaps.hpp generateindexmaps.cpp 
	$(CXX) generateindexmaps.cpp -c -o generateindexmaps.o 

generatemultiindices.o: generatemultiindices.hpp generatemultiindices.cpp 
	$(CXX) generatemultiindices.cpp -c -o generatemultiindices.o 

heappermgen.o: heappermgen.hpp heappermgen.cpp 
	$(CXX) heappermgen.cpp -c -o heappermgen.o 

objects = indexrange.o indexmap.o multiindex.o generateindexmaps.o generatemultiindices.o heappermgen.o

    
    

    
all: $(objects) 

include ../makefile.clean
