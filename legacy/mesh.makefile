SHELL = /bin/sh

default: all

objects = mesh.o coordinates.o io.coordinates.o mesh.simplicial1D.o io.simplicial1D.o mesh.simplicial2D.o io.simplicial2D.o mesh.simplicial3D.o io.simplicial3D.o mesh.simplicialND.o io.simplicialND.o

all: $(objects) libmesh.so 

include ../common.make 

include ../makefile.clean

libmesh.so: $(objects)
	g++ -shared -o libmesh.so $(objects)

coordinates.o: coordinates.cpp coordinates.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) coordinates.cpp -c -o coordinates.o 

io.coordinates.o: io.coordinates.cpp io.coordinates.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.coordinates.cpp -c -o io.coordinates.o 

mesh.o: mesh.cpp mesh.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.cpp -c -o mesh.o 

mesh.simplicial1D.o: mesh.simplicial1D.cpp mesh.simplicial1D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.simplicial1D.cpp -c -o mesh.simplicial1D.o 

io.simplicial1D.o: io.simplicial1D.cpp io.simplicial1D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.simplicial1D.cpp -c -o io.simplicial1D.o 

mesh.simplicial2D.o: mesh.simplicial2D.cpp mesh.simplicial2D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.simplicial2D.cpp -c -o mesh.simplicial2D.o 

io.simplicial2D.o: io.simplicial2D.cpp io.simplicial2D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.simplicial2D.cpp -c -o io.simplicial2D.o 

mesh.simplicial3D.o: mesh.simplicial3D.cpp mesh.simplicial3D.hpp mesh.simplicial3D.br.cpp mesh.simplicial3D.ur.cpp mesh.simplicial3D.mpr.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.simplicial3D.cpp -c -o mesh.simplicial3D.o 

io.simplicial3D.o: io.simplicial3D.cpp io.simplicial3D.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.simplicial3D.cpp -c -o io.simplicial3D.o 

mesh.simplicialND.o: mesh.simplicialND.cpp mesh.simplicialND.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.simplicialND.cpp -c -o mesh.simplicialND.o 

io.simplicialND.o: io.simplicialND.cpp io.simplicialND.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.simplicialND.cpp -c -o io.simplicialND.o 

# mesh.manifold2D.o: mesh.manifold2D.cpp mesh.manifold2D.hpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh.manifold2D.cpp -c -o mesh.manifold2D.o 
# 
# io.manifold2D.o: io.manifold2D.cpp io.manifold2D.hpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.manifold2D.cpp -c -o io.manifold2D.o 
# 
# vtkwriter.manifold2D.o: vtkwriter.manifold2D.cpp vtkwriter.manifold2D.hpp
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtkwriter.manifold2D.cpp -c -o vtkwriter.manifold2D.o  

# io.simplicialmesh.o simplicialmesh.o io.simplicialmesh.o generatesimplicialmesh.o vtkwriter.o manifold.2D.o 
# mesh.manifold2D.o io.manifold2D.o


