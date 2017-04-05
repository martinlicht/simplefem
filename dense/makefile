
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

dense.matrixtensorproduct.o: dense.matrixtensorproduct.cpp dense.matrixtensorproduct.hpp ../operators/floatvector.hpp
	$(CC) dense.matrixtensorproduct.cpp -c -o dense.matrixtensorproduct.o 

dense.functions.o: dense.functions.cpp dense.functions.hpp ../operators/floatvector.hpp
	$(CC) dense.functions.cpp -c -o dense.functions.o 

dense.scalarfunctions.o: dense.scalarfunctions.cpp dense.scalarfunctions.hpp ../operators/floatvector.hpp
	$(CC) dense.scalarfunctions.cpp -c -o dense.scalarfunctions.o 

dense.simplesolver.o: dense.simplesolver.cpp dense.simplesolver.hpp ../operators/floatvector.hpp
	$(CC) dense.simplesolver.cpp -c -o dense.simplesolver.o 

qr.factorization.o: qr.factorization.cpp qr.factorization.hpp ../operators/floatvector.hpp
	$(CC) qr.factorization.cpp -c -o qr.factorization.o 

lu.factorization.o: lu.factorization.cpp lu.factorization.hpp ../operators/floatvector.hpp
	$(CC) lu.factorization.cpp -c -o lu.factorization.o 

cholesky.o: cholesky.cpp cholesky.hpp ../operators/floatvector.hpp
	$(CC) cholesky.cpp -c -o cholesky.o 

objects = densematrix.o dense.functions.o dense.scalarfunctions.o dense.simplesolver.o qr.factorization.o lu.factorization.o cholesky.o dense.matrixtensorproduct.o

all: $(objects)

include ../makefile.clean
