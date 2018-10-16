SHELL = /bin/sh

default: all

objects = vtkwriter.mesh1D.o vtkwriter.mesh2D.o vtkwriter.mesh3D.o   

all: $(objects)

include ../makefile.rules 

include ../makefile.clean

vtkwriter.mesh1D.o: vtkwriter.mesh1D.hpp vtkwriter.mesh1D.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh1D.cpp -c -o vtkwriter.mesh1D.o  

vtkwriter.mesh2D.o: vtkwriter.mesh2D.hpp vtkwriter.mesh2D.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh2D.cpp -c -o vtkwriter.mesh2D.o  

vtkwriter.mesh3D.o: vtkwriter.mesh3D.hpp vtkwriter.mesh3D.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh3D.cpp -c -o vtkwriter.mesh3D.o  
