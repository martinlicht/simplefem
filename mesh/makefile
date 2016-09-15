
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp
	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

generatesimplicialmesh.o: generatesimplicialmesh.hpp generatesimplicialmesh.cpp
	$(CC) generatesimplicialmesh.cpp -c -o generatesimplicialmesh.o 

io.simplicialmesh.o: io.simplicialmesh.hpp io.simplicialmesh.cpp
	$(CC) io.simplicialmesh.cpp -c -o io.simplicialmesh.o 

io.coordinates.o: io.coordinates.hpp io.coordinates.cpp
	$(CC) io.coordinates.cpp -c -o io.coordinates.o 

vtkwriter.o: vtkwriter.hpp vtkwriter.cpp
	$(CC) vtkwriter.cpp -c -o vtkwriter.o  

objects = coordinates.o simplicialmesh.o io.simplicialmesh.o io.coordinates.o generatesimplicialmesh.o vtkwriter.o io.coordinates.o io.simplicialmesh.o 



unittest.coordinates.out: unittest.coordinates.cpp coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp coordinates.o ../combinatorics/*.o ../operators/*.o -o unittest.coordinates.out

unittest.simplicialmesh.out: unittest.simplicialmesh.cpp simplicialmesh.o simplicialmesh.hpp unittest.simplicialmesh.cpp generatesimplicialmesh.o vtkwriter.o
	$(CC) unittest.simplicialmesh.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o vtkwriter.o ../operators/*.o ../combinatorics/*.o -o unittest.simplicialmesh.out

unittest.vtkwriter.out: unittest.vtkwriter.cpp simplicialmesh.o generatesimplicialmesh.o vtkwriter.o
	$(CC) unittest.vtkwriter.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o vtkwriter.o ../operators/*.o ../combinatorics/*.o -o unittest.vtkwriter.out

unittest.io.out: unittest.io.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o io.coordinates.o io.simplicialmesh.o 
	$(CC) unittest.io.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o io.coordinates.o io.simplicialmesh.o ../operators/*.o ../combinatorics/*.o -o unittest.io.out

unittests = unittest.coordinates.out unittest.simplicialmesh.out unittest.io.out unittest.vtkwriter.out



all: $(objects) $(unittests) 


clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out *.temp