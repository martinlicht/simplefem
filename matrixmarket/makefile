SHELL = /bin/sh

default: all

objects = matrixmarket.o 

all: $(objects) libmatrixmarket.so

include ../common.make 

include ../makefile.clean

libmatrixmarket.so: $(objects)
	g++ -shared -o libmatrixmarket.so $(objects)

matrixmarket.o: matrixmarket.cpp matrixmarket.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matrixmarket.cpp -c -o matrixmarket.o 
