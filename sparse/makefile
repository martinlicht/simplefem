
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

sparsematrix.o: sparsematrix.hpp sparsematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) sparsematrix.cpp -c -o sparsematrix.o 

objects = sparsematrix.o 


unittest.sparsematrix.out: unittest.sparsematrix.cpp sparsematrix.o $(objects)
	$(CC) unittest.sparsematrix.cpp ../combinatorics/*.o *.o -o unittest.sparsematrix.out

unittests = unittest.sparsematrix.out 

all: $(objects) $(unittests)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 