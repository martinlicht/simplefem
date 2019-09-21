#!/usr/bin/make

SHELL = /bin/sh

default: all

objects = indexrange.o indexmap.o multiindex.o generateindexmaps.o generatemultiindices.o heappermgen.o

all: $(objects) libcombinatorics.so

include ../common.make 

include ../makefile.clean

libcombinatorics.so: $(objects)
	g++ -shared -o libcombinatorics.so $(objects)

indexrange.o: indexrange.cpp indexrange.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) indexrange.cpp -c -o indexrange.o 

indexmap.o: indexmap.cpp indexmap.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) indexmap.cpp -c -o indexmap.o 

multiindex.o: multiindex.cpp multiindex.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) multiindex.cpp -c -o multiindex.o 

generateindexmaps.o: generateindexmaps.cpp generateindexmaps.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) generateindexmaps.cpp -c -o generateindexmaps.o 

generatemultiindices.o: generatemultiindices.cpp generatemultiindices.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) generatemultiindices.cpp -c -o generatemultiindices.o 

heappermgen.o: heappermgen.cpp heappermgen.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) heappermgen.cpp -c -o heappermgen.o 
