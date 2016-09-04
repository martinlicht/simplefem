
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



massmatrix.element.o: massmatrix.element.hpp massmatrix.element.cpp 
	$(CC) massmatrix.element.cpp -c -o massmatrix.element.o 

objects = massmatrix.element.o 

    
    
unittest.massmatrix.element.out: massmatrix.element.o massmatrix.element.hpp unittest.massmatrix.element.cpp 
	$(CC) unittest.massmatrix.element.cpp massmatrix.element.o -o unittest.massmatrix.element.out

# unittests = unittest.massmatrix.element.out
unittests = 
    

    
all: $(objects) $(unittests) 

clean: 
	rm -f *.o *.exe *.exe.stackdump *.gch *.out