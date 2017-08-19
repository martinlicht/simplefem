
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



vtkwriter.mesh2D.o: vtkwriter.mesh2D.hpp vtkwriter.mesh2D.cpp
	$(CC) vtkwriter.mesh2D.cpp -c -o vtkwriter.mesh2D.o  

objects = vtkwriter.mesh2D.o   



all: $(objects) $(unittests) 


include ../makefile.clean
