
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



matrixmarket.o: matrixmarket.hpp matrixmarket.cpp 
	$(CC) matrixmarket.cpp -c -o matrixmarket.o 

objects = matrixmarket.o 

all: $(objects) 

include ../makefile.clean
