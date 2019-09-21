SHELL = /bin/sh

default: all

objects = vtkwriter.mesh1D.o vtkwriter.mesh2D.o vtkwriter.mesh3D.o   

all: $(objects) libvtk.so

include ../common.make 

include ../makefile.clean

libvtk.so: $(objects)
	g++ -shared -o libvtk.so $(objects)

vtkwriter.mesh1D.o: vtkwriter.mesh1D.cpp vtkwriter.mesh1D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh1D.cpp -c -o vtkwriter.mesh1D.o  

vtkwriter.mesh2D.o: vtkwriter.mesh2D.cpp vtkwriter.mesh2D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh2D.cpp -c -o vtkwriter.mesh2D.o  

vtkwriter.mesh3D.o: vtkwriter.mesh3D.cpp vtkwriter.mesh3D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.mesh3D.cpp -c -o vtkwriter.mesh3D.o  
