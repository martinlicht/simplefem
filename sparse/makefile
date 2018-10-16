SHELL = /bin/sh

default: all

objects = sparsematrix.o readwrite.o 

all: $(objects) sparse.so

include ../makefile.rules 

include ../makefile.clean

sparse.so: $(objects)
	g++ -shared -o sparse.so $(objects)

sparsematrix.o: sparsematrix.hpp sparsematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) sparsematrix.cpp -c -o sparsematrix.o 

readwrite.o: readwrite.hpp readwrite.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) readwrite.cpp -c -o readwrite.o 
