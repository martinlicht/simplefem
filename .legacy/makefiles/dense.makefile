SHELL = /bin/sh

default: all

objects = densematrix.o readwrite.o functions.o scalarfunctions.o simplesolver.o qr.factorization.o gaussjordan.o cholesky.o matrixtensorproduct.o iterativeinverse.o

all: $(objects) libdense.so

include ../common.make 

include ../makefile.clean

libdense.so: $(objects)
	g++ -shared -o libdense.so $(objects)

densematrix.o: densematrix.cpp densematrix.hpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) densematrix.cpp -c -o densematrix.o 

readwrite.o: readwrite.cpp readwrite.hpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
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

gaussjordan.o: gaussjordan.cpp gaussjordan.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) gaussjordan.cpp -c -o gaussjordan.o 

cholesky.o: cholesky.cpp cholesky.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) cholesky.cpp -c -o cholesky.o 

iterativeinverse.o: iterativeinverse.cpp iterativeinverse.hpp ../operators/floatvector.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) iterativeinverse.cpp -c -o iterativeinverse.o 
