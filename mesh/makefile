
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



coordinates.o: coordinates.hpp coordinates.cpp 
	$(CC) coordinates.cpp -c -o coordinates.o 

simplicialmesh.o: simplicialmesh.hpp simplicialmesh.cpp coordinates.o
	$(CC) simplicialmesh.cpp -c -o simplicialmesh.o 

generatesimplicialmesh.o: generatesimplicialmesh.hpp generatesimplicialmesh.cpp simplicialmesh.o coordinates.o
	$(CC) generatesimplicialmesh.cpp -c -o generatesimplicialmesh.o 

io.simplicialmesh.o: io.simplicialmesh.hpp io.simplicialmesh.cpp simplicialmesh.o coordinates.o
	$(CC) io.simplicialmesh.cpp -c -o io.simplicialmesh.o 

io.coordinates.o: io.coordinates.hpp io.coordinates.cpp coordinates.o simplicialmesh.o
	$(CC) io.coordinates.cpp -c -o io.coordinates.o 

vtkwriter.o: vtkwriter.hpp vtkwriter.cpp coordinates.o simplicialmesh.o
	$(CC) vtkwriter.cpp -c -o vtkwriter.o  

objects = coordinates.o simplicialmesh.o io.simplicialmesh.o io.coordinates.o generatesimplicialmesh.o vtkwriter.o 



unittest.coordinates.out: coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp coordinates.o ../combinatorics/*.o ../operators/*.o -o unittest.coordinates.out

unittest.simplicialmesh.out: simplicialmesh.o simplicialmesh.hpp unittest.simplicialmesh.cpp generatesimplicialmesh.o vtkwriter.o
	$(CC) unittest.simplicialmesh.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o vtkwriter.o ../operators/*.o ../combinatorics/*.o -o unittest.simplicialmesh.out

unittests = unittest.coordinates.out unittest.simplicialmesh.out 



all: $(objects) $(unittests) 


clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 