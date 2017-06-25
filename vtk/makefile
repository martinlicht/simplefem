
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



vtkwriter.manifold2D.o: vtkwriter.manifold2D.hpp vtkwriter.manifold2D.cpp
	$(CC) vtkwriter.manifold2D.cpp -c -o vtkwriter.manifold2D.o  

objects = vtkwriter.manifold2D.o   



all: $(objects) $(unittests) 


include ../makefile.clean
