
# CXX = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) densematrix.cpp -c -o densematrix.o 

readwrite.o: readwrite.hpp readwrite.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) readwrite.cpp -c -o readwrite.o 

matrixtensorproduct.o: matrixtensorproduct.cpp matrixtensorproduct.hpp ../operators/floatvector.hpp
	$(CXX) matrixtensorproduct.cpp -c -o matrixtensorproduct.o 

functions.o: functions.cpp functions.hpp ../operators/floatvector.hpp
	$(CXX) functions.cpp -c -o functions.o 

scalarfunctions.o: scalarfunctions.cpp scalarfunctions.hpp ../operators/floatvector.hpp
	$(CXX) scalarfunctions.cpp -c -o scalarfunctions.o 

simplesolver.o: simplesolver.cpp simplesolver.hpp ../operators/floatvector.hpp
	$(CXX) simplesolver.cpp -c -o simplesolver.o 

qr.factorization.o: qr.factorization.cpp qr.factorization.hpp ../operators/floatvector.hpp
	$(CXX) qr.factorization.cpp -c -o qr.factorization.o 

lu.factorization.o: lu.factorization.cpp lu.factorization.hpp ../operators/floatvector.hpp
	$(CXX) lu.factorization.cpp -c -o lu.factorization.o 

cholesky.o: cholesky.cpp cholesky.hpp ../operators/floatvector.hpp
	$(CXX) cholesky.cpp -c -o cholesky.o 

iterativeinverse.o: iterativeinverse.cpp iterativeinverse.hpp ../operators/floatvector.hpp
	$(CXX) iterativeinverse.cpp -c -o iterativeinverse.o 

objects = densematrix.o readwrite.o functions.o scalarfunctions.o simplesolver.o qr.factorization.o lu.factorization.o cholesky.o matrixtensorproduct.o iterativeinverse.o

all: $(objects)

include ../makefile.clean
