
# CXX = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



massmatrix.element.o: massmatrix.element.hpp massmatrix.element.cpp 
	$(CXX) massmatrix.element.cpp -c -o massmatrix.element.o 

diffmatrix.element.o: diffmatrix.element.hpp diffmatrix.element.cpp 
	$(CXX) diffmatrix.element.cpp -c -o diffmatrix.element.o 

objects = massmatrix.element.o diffmatrix.element.o

    
    
unittest.massmatrix.element.out: massmatrix.element.o massmatrix.element.hpp unittest.massmatrix.element.cpp 
	$(CXX) unittest.massmatrix.element.cpp massmatrix.element.o -o unittest.massmatrix.element.out

# unittests = unittest.massmatrix.element.out
unittests = 
    

    
all: $(objects) $(unittests) 

include ../makefile.clean
