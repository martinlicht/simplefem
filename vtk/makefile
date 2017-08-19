
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



vtkwriter.mesh1D.o: vtkwriter.mesh1D.hpp vtkwriter.mesh1D.cpp
	$(CC) vtkwriter.mesh1D.cpp -c -o vtkwriter.mesh1D.o  

vtkwriter.mesh2D.o: vtkwriter.mesh2D.hpp vtkwriter.mesh2D.cpp
	$(CC) vtkwriter.mesh2D.cpp -c -o vtkwriter.mesh2D.o  

vtkwriter.mesh3D.o: vtkwriter.mesh3D.hpp vtkwriter.mesh3D.cpp
	$(CC) vtkwriter.mesh3D.cpp -c -o vtkwriter.mesh3D.o  

objects = vtkwriter.mesh1D.o vtkwriter.mesh2D.o vtkwriter.mesh3D.o   



all: $(objects) $(unittests) 


include ../makefile.clean
