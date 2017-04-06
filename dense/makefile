
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

matrixtensorproduct.o: matrixtensorproduct.cpp matrixtensorproduct.hpp ../operators/floatvector.hpp
	$(CC) matrixtensorproduct.cpp -c -o matrixtensorproduct.o 

functions.o: functions.cpp functions.hpp ../operators/floatvector.hpp
	$(CC) functions.cpp -c -o functions.o 

scalarfunctions.o: scalarfunctions.cpp scalarfunctions.hpp ../operators/floatvector.hpp
	$(CC) scalarfunctions.cpp -c -o scalarfunctions.o 

simplesolver.o: simplesolver.cpp simplesolver.hpp ../operators/floatvector.hpp
	$(CC) simplesolver.cpp -c -o simplesolver.o 

qr.factorization.o: qr.factorization.cpp qr.factorization.hpp ../operators/floatvector.hpp
	$(CC) qr.factorization.cpp -c -o qr.factorization.o 

lu.factorization.o: lu.factorization.cpp lu.factorization.hpp ../operators/floatvector.hpp
	$(CC) lu.factorization.cpp -c -o lu.factorization.o 

cholesky.o: cholesky.cpp cholesky.hpp ../operators/floatvector.hpp
	$(CC) cholesky.cpp -c -o cholesky.o 

iterativeinverse.o: iterativeinverse.cpp iterativeinverse.hpp ../operators/floatvector.hpp
	$(CC) iterativeinverse.cpp -c -o iterativeinverse.o 

objects = densematrix.o functions.o scalarfunctions.o simplesolver.o qr.factorization.o lu.factorization.o cholesky.o matrixtensorproduct.o iterativeinverse.o

all: $(objects)

include ../makefile.clean
