SHELL = /bin/sh

default: all

objects = sparsematrix.o matcsr.o readwrite.o 

all: $(objects) libsparse.so

include ../common.make 

include ../makefile.clean

libsparse.so: $(objects)
	g++ -shared -o libsparse.so $(objects)

matcsr.o: matcsr.cpp matcsr.hpp sparsematrix.hpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matcsr.cpp -c -o matcsr.o 

sparsematrix.o: sparsematrix.cpp sparsematrix.hpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) sparsematrix.cpp -c -o sparsematrix.o 

readwrite.o: readwrite.cpp readwrite.hpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) readwrite.cpp -c -o readwrite.o 
