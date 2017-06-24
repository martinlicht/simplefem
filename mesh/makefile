
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

# # simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp
# # 	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

# # manifold.2D.o: manifold.2D.hpp manifold.2D.cpp
# # 	$(CC) manifold.2D.cpp -c -o manifold.2D.o 

mesh.o: mesh.hpp mesh.cpp
	$(CC) mesh.cpp -c -o mesh.o 

mesh.manifold2D.o: mesh.manifold2D.hpp mesh.manifold2D.cpp
	$(CC) mesh.manifold2D.cpp -c -o mesh.manifold2D.o 

# generatesimplicialmesh.o: generatesimplicialmesh.hpp generatesimplicialmesh.cpp
# 	$(CC) generatesimplicialmesh.cpp -c -o generatesimplicialmesh.o 

# io.simplicialmesh.o: io.simplicialmesh.hpp io.simplicialmesh.cpp
# 	$(CC) io.simplicialmesh.cpp -c -o io.simplicialmesh.o 

io.coordinates.o: io.coordinates.hpp io.coordinates.cpp
	$(CC) io.coordinates.cpp -c -o io.coordinates.o 

# vtkwriter.o: vtkwriter.hpp vtkwriter.cpp
# 	$(CC) vtkwriter.cpp -c -o vtkwriter.o  

vtkwriter.manifold2D.o: vtkwriter.manifold2D.hpp vtkwriter.manifold2D.cpp
	$(CC) vtkwriter.manifold2D.cpp -c -o vtkwriter.manifold2D.o  

objects = coordinates.o mesh.o mesh.manifold2D.o io.coordinates.o vtkwriter.manifold2D.o io.coordinates.o  
# io.simplicialmesh.o simplicialmesh.o io.simplicialmesh.o generatesimplicialmesh.o vtkwriter.o manifold.2D.o 



all: $(objects) $(unittests) 


include ../makefile.clean
