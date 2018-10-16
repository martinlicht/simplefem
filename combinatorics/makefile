SHELL = /bin/sh

default: all

objects = indexrange.o indexmap.o multiindex.o generateindexmaps.o generatemultiindices.o heappermgen.o

all: $(objects) combinatorics.so

include ../makefile.rules 

include ../makefile.clean

combinatorics.so: $(objects)
	g++ -shared -o combinatorics.so $(objects)

indexrange.o: indexrange.hpp indexrange.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.hpp indexmap.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) indexmap.cpp -c -o indexmap.o 

multiindex.o: multiindex.hpp multiindex.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) multiindex.cpp -c -o multiindex.o 

generateindexmaps.o: generateindexmaps.hpp generateindexmaps.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) generateindexmaps.cpp -c -o generateindexmaps.o 

generatemultiindices.o: generatemultiindices.hpp generatemultiindices.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) generatemultiindices.cpp -c -o generatemultiindices.o 

heappermgen.o: heappermgen.hpp heappermgen.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) heappermgen.cpp -c -o heappermgen.o 
