# TODO: This makefile still needs to be written properly

SHELL = /bin/sh

default: all

include ../makefile.rules 



massmatrix.element.o: massmatrix.element.hpp massmatrix.element.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) massmatrix.element.cpp -c -o massmatrix.element.o 

diffmatrix.element.o: diffmatrix.element.hpp diffmatrix.element.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) diffmatrix.element.cpp -c -o diffmatrix.element.o 

objects = massmatrix.element.o diffmatrix.element.o

    
    
unittest.massmatrix.element.out: massmatrix.element.o massmatrix.element.hpp unittest.massmatrix.element.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) unittest.massmatrix.element.cpp massmatrix.element.o -o unittest.massmatrix.element.out

# unittests = unittest.massmatrix.element.out
unittests = 
    

    
all: $(objects) $(unittests) 

include ../makefile.clean
