SHELL = /bin/sh

default: all

objects = matrixmarket.o 

all: $(objects) 

include ../makefile.rules 

include ../makefile.clean

matrixmarket.o: matrixmarket.hpp matrixmarket.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matrixmarket.cpp -c -o matrixmarket.o 
