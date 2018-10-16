SHELL = /bin/sh

default: all

objects = matrixmarket.o 

all: $(objects) matrixmarket.so

include ../makefile.rules 

include ../makefile.clean

matrixmarket.so: $(objects)
	g++ -shared -o matrixmarket.so $(objects)

matrixmarket.o: matrixmarket.hpp matrixmarket.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matrixmarket.cpp -c -o matrixmarket.o 
