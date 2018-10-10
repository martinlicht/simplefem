
# CXX = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) densematrix.cpp -c -o densematrix.o 

readwrite.o: readwrite.hpp readwrite.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) readwrite.cpp -c -o readwrite.o 

matrixtensorproduct.o: matrixtensorproduct.cpp matrixtensorproduct.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matrixtensorproduct.cpp -c -o matrixtensorproduct.o 

functions.o: functions.cpp functions.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) functions.cpp -c -o functions.o 

scalarfunctions.o: scalarfunctions.cpp scalarfunctions.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) scalarfunctions.cpp -c -o scalarfunctions.o 

simplesolver.o: simplesolver.cpp simplesolver.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) simplesolver.cpp -c -o simplesolver.o 

qr.factorization.o: qr.factorization.cpp qr.factorization.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) qr.factorization.cpp -c -o qr.factorization.o 

lu.factorization.o: lu.factorization.cpp lu.factorization.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) lu.factorization.cpp -c -o lu.factorization.o 

cholesky.o: cholesky.cpp cholesky.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) cholesky.cpp -c -o cholesky.o 

iterativeinverse.o: iterativeinverse.cpp iterativeinverse.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) iterativeinverse.cpp -c -o iterativeinverse.o 

objects = densematrix.o readwrite.o functions.o scalarfunctions.o simplesolver.o qr.factorization.o lu.factorization.o cholesky.o matrixtensorproduct.o iterativeinverse.o

all: $(objects)

include ../makefile.clean
