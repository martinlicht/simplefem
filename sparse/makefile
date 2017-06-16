
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

sparsematrix.o: sparsematrix.hpp sparsematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) sparsematrix.cpp -c -o sparsematrix.o 

readwrite.o: readwrite.hpp readwrite.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) readwrite.cpp -c -o readwrite.o 

objects = sparsematrix.o readwrite.o 



all: $(objects) 

include ../makefile.clean
