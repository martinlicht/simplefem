
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all



unittest.indexrange.out: indexrange.o indexrange.hpp unittest.indexrange.cpp 
	$(CC) unittest.indexrange.cpp indexrange.o -o unittest.indexrange.out

unittest.indexmap.out: indexmap.o indexrange.o indexmap.hpp unittest.indexmap.cpp 
	$(CC) unittest.indexmap.cpp indexmap.o indexrange.o -o unittest.indexmap.out

unittest.multiindex.out: multiindex.o indexrange.o indexmap.hpp unittest.multiindex.cpp 
	$(CC) unittest.multiindex.cpp indexmap.o indexrange.o multiindex.o -o unittest.multiindex.out

unittest.generateindexmaps.out: generateindexmaps.o indexmap.o indexrange.o indexmap.hpp unittest.generateindexmaps.cpp 
	$(CC) unittest.generateindexmaps.cpp generateindexmaps.o indexmap.o indexrange.o -o unittest.generateindexmaps.out

unittest.generatemultiindices.out: generatemultiindices.o indexmap.o indexrange.o multiindex.o unittest.generatemultiindices.cpp 
	$(CC) unittest.generatemultiindices.cpp generatemultiindices.o multiindex.o indexmap.o indexrange.o -o unittest.generatemultiindices.out

combinatorics.unittests = unittest.indexrange.out unittest.indexmap.out unittest.multiindex.out unittest.generateindexmaps.out unittest.generatemultiindices.out
    



unittest.floatvector.out: floatvector.o floatvector.hpp unittest.floatvector.cpp $(objects) 
	$(CC) unittest.floatvector.cpp floatvector.o -o unittest.floatvector.out

unittest.scalingoperator.out: unittest.scalingoperator.cpp scalingoperator.o $(objects)
	$(CC) unittest.scalingoperator.cpp scalingoperator.o floatvector.o linearoperator.o -o unittest.scalingoperator.out
	
unittest.diagonaloperator.out: unittest.diagonaloperator.cpp diagonaloperator.o $(objects)
	$(CC) unittest.diagonaloperator.cpp diagonaloperator.o scalingoperator.o floatvector.o linearoperator.o -o unittest.diagonaloperator.out
	
operators.unittests = unittest.floatvector.out unittest.scalingoperator.out unittest.diagonaloperator.out 


    
    
    
    
unittest.densematrix.out: unittest.densematrix.cpp densematrix.o $(objects) 
	$(CC) unittest.densematrix.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o unittest.densematrix.out 

unittest.matrixalgorithm.out: unittest.matrixalgorithm.cpp matrixalgorithm.o $(objects) 
	$(CC) unittest.matrixalgorithm.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o unittest.matrixalgorithm.out

dense.unittests = unittest.densematrix.out unittest.matrixalgorithm.out




unittest.sparsematrix.out: unittest.sparsematrix.cpp sparsematrix.o $(objects)
	$(CC) unittest.sparsematrix.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o *.o -o unittest.sparsematrix.out

sparse.unittests = unittest.sparsematrix.out 



unittest.crm.out: unittest.crm.cpp crm.o $(objects)
	$(CC) unittest.crm.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o unittest.crm.out

moretesting.out: moretesting.cpp crm.o $(objects)
	$(CC) moretesting.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o moretesting.out

solver.unittests = unittest.crm.out moretesting.out 





unittest.coordinates.out: unittest.coordinates.cpp coordinates.o coordinates.hpp unittest.coordinates.cpp 
	$(CC) unittest.coordinates.cpp coordinates.o ../combinatorics/*.o ../operators/*.o -o unittest.coordinates.out

unittest.simplicialmesh.out: unittest.simplicialmesh.cpp simplicialmesh.o simplicialmesh.hpp unittest.simplicialmesh.cpp generatesimplicialmesh.o vtkwriter.o
	$(CC) unittest.simplicialmesh.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o vtkwriter.o ../operators/*.o ../combinatorics/*.o -o unittest.simplicialmesh.out

unittest.vtkwriter.out: unittest.vtkwriter.cpp simplicialmesh.o generatesimplicialmesh.o vtkwriter.o
	$(CC) unittest.vtkwriter.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o vtkwriter.o ../operators/*.o ../combinatorics/*.o -o unittest.vtkwriter.out

unittest.io.out: unittest.io.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o io.coordinates.o io.simplicialmesh.o 
	$(CC) unittest.io.cpp coordinates.o simplicialmesh.o generatesimplicialmesh.o io.coordinates.o io.simplicialmesh.o ../operators/*.o ../combinatorics/*.o -o unittest.io.out

mesh.unittests = unittest.coordinates.out unittest.simplicialmesh.out unittest.io.out unittest.vtkwriter.out




    
all: $(combinatorics.unittests) $(operators.unittests) $(dense.unittests) $(sparse.unittests) $(solver.unittests) $(mesh.unittests)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 