
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

Objects = ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o ../solver/*.o ../mesh/*.o
Headers = ../combinatorics/*.hpp ../operators/*.hpp ../dense/*.hpp ../sparse/*.hpp ../solver/*.hpp ../mesh/*.hpp

Includes = $(Headers) $(Objects) 

test.combinatorics.indexrange.out: test.combinatorics.indexrange.cpp $(Includes) 
	$(CC) test.combinatorics.indexrange.cpp ../combinatorics/*.o -o test.combinatorics.indexrange.out

test.combinatorics.indexmap.out: test.combinatorics.indexmap.cpp $(Includes) 
	$(CC) test.combinatorics.indexmap.cpp ../combinatorics/*.o -o test.combinatorics.indexmap.out

test.combinatorics.multiindex.out: test.combinatorics.multiindex.cpp $(Includes) 
	$(CC) test.combinatorics.multiindex.cpp ../combinatorics/*.o -o test.combinatorics.multiindex.out

test.combinatorics.generateindexmaps.out: test.combinatorics.generateindexmaps.cpp $(Includes) 
	$(CC) test.combinatorics.generateindexmaps.cpp ../combinatorics/*.o -o test.combinatorics.generateindexmaps.out

test.combinatorics.generatemultiindices.out: test.combinatorics.generatemultiindices.cpp $(Includes) 
	$(CC) test.combinatorics.generatemultiindices.cpp ../combinatorics/*.o -o test.combinatorics.generatemultiindices.out

combinatorics.tests = test.combinatorics.indexrange.out test.combinatorics.indexmap.out test.combinatorics.multiindex.out test.combinatorics.generateindexmaps.out test.combinatorics.generatemultiindices.out
    



test.operators.floatvector.out: test.operators.floatvector.cpp $(Includes) 
	$(CC) test.operators.floatvector.cpp $(Objects) -o test.operators.floatvector.out

test.operators.scalingoperator.out: test.operators.scalingoperator.cpp $(Includes)
	$(CC) test.operators.scalingoperator.cpp $(Objects) -o test.operators.scalingoperator.out
	
test.operators.diagonaloperator.out: test.operators.diagonaloperator.cpp $(Includes)
	$(CC) test.operators.diagonaloperator.cpp $(Objects) -o test.operators.diagonaloperator.out
	
test.operators.sumoperator.out: test.operators.sumoperator.cpp $(Includes)
	$(CC) test.operators.sumoperator.cpp $(Objects) -o test.operators.sumoperator.out
	
test.operators.productoperator.out: test.operators.productoperator.cpp $(Includes)
	$(CC) test.operators.productoperator.cpp $(Objects) -o test.operators.productoperator.out
	
operators.tests = test.operators.floatvector.out test.operators.scalingoperator.out test.operators.diagonaloperator.out test.operators.sumoperator.out test.operators.productoperator.out 


    
    
    
    
test.dense.densematrix.out: test.dense.densematrix.cpp $(Includes) 
	$(CC) test.dense.densematrix.cpp $(Objects) -o test.dense.densematrix.out 

test.dense.matrixalgorithm.out: test.dense.matrixalgorithm.cpp $(Includes) 
	$(CC) test.dense.matrixalgorithm.cpp $(Objects) -o test.dense.matrixalgorithm.out

dense.tests = test.dense.densematrix.out test.dense.matrixalgorithm.out




test.sparse.sparsematrix.out: test.sparse.sparsematrix.cpp $(Includes)
	$(CC) test.sparse.sparsematrix.cpp $(Objects) -o test.sparse.sparsematrix.out

sparse.tests = test.sparse.sparsematrix.out 



test.solver.crm.out: test.solver.crm.cpp $(Includes)
	$(CC) test.solver.crm.cpp $(Objects) -o test.solver.crm.out

test.solver.moretesting.out: test.solver.moretesting.cpp $(Includes)
	$(CC) test.solver.moretesting.cpp $(Objects) -o test.solver.moretesting.out

solver.tests = test.solver.crm.out test.solver.moretesting.out 





test.mesh.coordinates.out: test.mesh.coordinates.cpp $(Includes) 
	$(CC) test.mesh.coordinates.cpp $(Objects) -o test.mesh.coordinates.out

test.mesh.simplicialmesh.out: test.mesh.simplicialmesh.cpp $(Includes)
	$(CC) test.mesh.simplicialmesh.cpp $(Objects) -o test.mesh.simplicialmesh.out

test.mesh.vtkwriter.out: test.mesh.vtkwriter.cpp $(Includes)
	$(CC) test.mesh.vtkwriter.cpp $(Objects) -o test.mesh.vtkwriter.out

test.mesh.io.out: test.mesh.io.cpp $(Includes)
	$(CC) test.mesh.io.cpp $(Objects) -o test.mesh.io.out

mesh.tests = test.mesh.coordinates.out test.mesh.simplicialmesh.out test.mesh.io.out test.mesh.vtkwriter.out




    
all: $(combinatorics.tests) $(operators.tests) 

# $(dense.tests) $(sparse.tests) $(solver.tests) $(mesh.tests)

include ../makefile.clean


