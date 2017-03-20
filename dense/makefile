
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

dense.functions.o: dense.functions.cpp dense.functions.hpp ../operators/floatvector.hpp
	$(CC) dense.functions.cpp -c -o dense.functions.o 

dense.simplesolver.o: dense.simplesolver.cpp dense.simplesolver.hpp ../operators/floatvector.hpp
	$(CC) dense.simplesolver.cpp -c -o dense.simplesolver.o 

dense.factorization.o: dense.factorization.cpp densematrix.hpp ../operators/floatvector.hpp
	$(CC) dense.factorization.cpp -c -o dense.factorization.o 

objects = densematrix.o dense.functions.o dense.simplesolver.o dense.factorization.o

all: $(objects)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 