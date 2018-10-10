
# CXX = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



coordinates.o: coordinates.hpp coordinates.cpp 
	$(CXX) coordinates.cpp -c -o coordinates.o 

io.coordinates.o: io.coordinates.hpp io.coordinates.cpp
	$(CXX) io.coordinates.cpp -c -o io.coordinates.o 

mesh.o: mesh.hpp mesh.cpp
	$(CXX) mesh.cpp -c -o mesh.o 

mesh.simplicial1D.o: mesh.simplicial1D.hpp mesh.simplicial1D.cpp
	$(CXX) mesh.simplicial1D.cpp -c -o mesh.simplicial1D.o 

io.simplicial1D.o: io.simplicial1D.hpp io.simplicial1D.cpp
	$(CXX) io.simplicial1D.cpp -c -o io.simplicial1D.o 

mesh.simplicial2D.o: mesh.simplicial2D.hpp mesh.simplicial2D.cpp
	$(CXX) mesh.simplicial2D.cpp -c -o mesh.simplicial2D.o 

io.simplicial2D.o: io.simplicial2D.hpp io.simplicial2D.cpp
	$(CXX) io.simplicial2D.cpp -c -o io.simplicial2D.o 

mesh.simplicial3D.o: mesh.simplicial3D.hpp mesh.simplicial3D.cpp
	$(CXX) mesh.simplicial3D.cpp -c -o mesh.simplicial3D.o 

io.simplicial3D.o: io.simplicial3D.hpp io.simplicial3D.cpp
	$(CXX) io.simplicial3D.cpp -c -o io.simplicial3D.o 

mesh.simplicialND.o: mesh.simplicialND.hpp mesh.simplicialND.cpp
	$(CXX) mesh.simplicialND.cpp -c -o mesh.simplicialND.o 

io.simplicialND.o: io.simplicialND.hpp io.simplicialND.cpp
	$(CXX) io.simplicialND.cpp -c -o io.simplicialND.o 

# mesh.manifold2D.o: mesh.manifold2D.hpp mesh.manifold2D.cpp
# 	$(CXX) mesh.manifold2D.cpp -c -o mesh.manifold2D.o 
# 
# io.manifold2D.o: io.manifold2D.hpp io.manifold2D.cpp
# 	$(CXX) io.manifold2D.cpp -c -o io.manifold2D.o 
# 
# vtkwriter.manifold2D.o: vtkwriter.manifold2D.hpp vtkwriter.manifold2D.cpp
# 	$(CXX) vtkwriter.manifold2D.cpp -c -o vtkwriter.manifold2D.o  

objects = mesh.o coordinates.o io.coordinates.o mesh.simplicial1D.o io.simplicial1D.o mesh.simplicial2D.o io.simplicial2D.o mesh.simplicialND.o io.simplicialND.o
# io.simplicialmesh.o simplicialmesh.o io.simplicialmesh.o generatesimplicialmesh.o vtkwriter.o manifold.2D.o  mesh.simplicial3D.o io.simplicial3D.o
# mesh.manifold2D.o io.manifold2D.o


all: $(objects) $(unittests) 


include ../makefile.clean
