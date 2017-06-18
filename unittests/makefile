
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

Objects.combinatorics = ../combinatorics/*.o
Objects.operators = ../operators/*.o
Objects.dense = ../dense/*.o
Objects.sparse = ../sparse/*.o
Objects.solver = ../solver/*.o
Objects.mesh = ../mesh/*.o

Headers.combinatorics = ../combinatorics/*.o
Headers.operators = ../operators/*.o
Headers.dense = ../dense/*.o
Headers.sparse = ../sparse/*.o
Headers.solver = ../solver/*.o
Headers.mesh = ../mesh/*.o

Includes.combinatorics = $(Headers.combinatorics) $(Objects.combinatorics) 
Includes.operators = $(Headers.operators) $(Objects.operators) $(Includes.combinatorics)
Includes.dense = $(Headers.dense) $(Objects.dense) $(Includes.operators)
Includes.sparse = $(Headers.sparse) $(Objects.sparse) $(Includes.operators)
Includes.solver = $(Headers.solver) $(Objects.solver) $(Includes.operators)
Includes.mesh = $(Headers.mesh) $(Objects.mesh) $(Objects.combinatorics) $(Includes.operators)


Objects = ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o ../solver/*.o ../mesh/*.o
Headers = ../combinatorics/*.hpp ../operators/*.hpp ../dense/*.hpp ../sparse/*.hpp ../solver/*.hpp ../mesh/*.hpp

Includes = $(Headers) $(Objects) 

test.combinatorics.indexrange.out: test.combinatorics.indexrange.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.indexrange.cpp ../combinatorics/*.o -o test.combinatorics.indexrange.out

test.combinatorics.indexmap.out: test.combinatorics.indexmap.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.indexmap.cpp ../combinatorics/*.o -o test.combinatorics.indexmap.out

test.combinatorics.multiindex.out: test.combinatorics.multiindex.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.multiindex.cpp ../combinatorics/*.o -o test.combinatorics.multiindex.out

test.combinatorics.generateindexmaps.out: test.combinatorics.generateindexmaps.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.generateindexmaps.cpp ../combinatorics/*.o -o test.combinatorics.generateindexmaps.out

test.combinatorics.generatemultiindices.out: test.combinatorics.generatemultiindices.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.generatemultiindices.cpp ../combinatorics/*.o -o test.combinatorics.generatemultiindices.out

test.combinatorics.heapsalgorithm.out: test.combinatorics.heapsalgorithm.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.heapsalgorithm.cpp ../combinatorics/*.o -o test.combinatorics.heapsalgorithm.out

combinatorics.tests = test.combinatorics.indexrange.out test.combinatorics.indexmap.out test.combinatorics.multiindex.out test.combinatorics.generateindexmaps.out test.combinatorics.generatemultiindices.out test.combinatorics.heapsalgorithm.out 
    



test.operators.floatvector.out: test.operators.floatvector.cpp $(Includes.operators) 
	$(CC) test.operators.floatvector.cpp $(Objects) -o test.operators.floatvector.out

test.operators.scalingoperator.out: test.operators.scalingoperator.cpp $(Includes.operators)
	$(CC) test.operators.scalingoperator.cpp $(Objects) -o test.operators.scalingoperator.out
	
test.operators.diagonaloperator.out: test.operators.diagonaloperator.cpp $(Includes.operators)
	$(CC) test.operators.diagonaloperator.cpp $(Objects) -o test.operators.diagonaloperator.out
	
test.operators.sumoperator.out: test.operators.sumoperator.cpp $(Includes.operators)
	$(CC) test.operators.sumoperator.cpp $(Objects) -o test.operators.sumoperator.out
	
test.operators.productoperator.out: test.operators.productoperator.cpp $(Includes.operators)
	$(CC) test.operators.productoperator.cpp $(Objects) -o test.operators.productoperator.out
	
test.operators.blockdiagonaloperator.out: test.operators.blockdiagonaloperator.cpp $(Includes.operators)
	$(CC) test.operators.blockdiagonaloperator.cpp $(Objects) -o test.operators.blockdiagonaloperator.out
	
test.operators.blockoperator.out: test.operators.blockoperator.cpp $(Includes.operators)
	$(CC) test.operators.blockoperator.cpp $(Objects) -o test.operators.blockoperator.out
	
operators.tests = test.operators.floatvector.out test.operators.scalingoperator.out test.operators.diagonaloperator.out test.operators.sumoperator.out test.operators.productoperator.out test.operators.blockdiagonaloperator.out test.operators.blockoperator.out


    
    
    
    
test.dense.densematrix.out: test.dense.densematrix.cpp $(Includes.dense) 
	$(CC) test.dense.densematrix.cpp $(Objects) -o test.dense.densematrix.out 

test.dense.transpose.out: test.dense.transpose.cpp $(Includes.dense) 
	$(CC) test.dense.transpose.cpp $(Objects) -o test.dense.transpose.out

test.dense.scalarfunctions.out: test.dense.scalarfunctions.cpp $(Includes.dense) 
	$(CC) test.dense.scalarfunctions.cpp $(Objects) -o test.dense.scalarfunctions.out

test.dense.detinvcof.out: test.dense.detinvcof.cpp $(Includes.dense) 
	$(CC) test.dense.detinvcof.cpp $(Objects) -o test.dense.detinvcof.out

test.dense.matrixalgorithm.out: test.dense.matrixalgorithm.cpp $(Includes.dense) 
	$(CC) test.dense.matrixalgorithm.cpp $(Objects) -o test.dense.matrixalgorithm.out

test.dense.cholesky.out: test.dense.cholesky.cpp $(Includes.dense) 
	$(CC) test.dense.cholesky.cpp $(Objects) -o test.dense.cholesky.out

test.dense.iterativeinverse.out: test.dense.iterativeinverse.cpp $(Includes.dense) 
	$(CC) test.dense.iterativeinverse.cpp $(Objects) -o test.dense.iterativeinverse.out

dense.tests = test.dense.densematrix.out test.dense.transpose.out test.dense.detinvcof.out test.dense.matrixalgorithm.out test.dense.cholesky.out test.dense.scalarfunctions.out test.dense.iterativeinverse.out




test.sparse.sparsematrix.out: test.sparse.sparsematrix.cpp $(Includes.sparse)
	$(CC) test.sparse.sparsematrix.cpp $(Objects) -o test.sparse.sparsematrix.out

sparse.tests = test.sparse.sparsematrix.out 



test.solver.crm.out: test.solver.crm.cpp $(Includes.solver)
	$(CC) test.solver.crm.cpp $(Objects) -o test.solver.crm.out

test.solver.moretesting.out: test.solver.moretesting.cpp $(Includes.solver)
	$(CC) test.solver.moretesting.cpp $(Objects) -o test.solver.moretesting.out

solver.tests = test.solver.crm.out test.solver.moretesting.out 





test.mesh.coordinates.out: test.mesh.coordinates.cpp $(Includes.mesh) 
	$(CC) test.mesh.coordinates.cpp $(Objects) -o test.mesh.coordinates.out

test.mesh.simplicialmesh.out: test.mesh.simplicialmesh.cpp $(Includes.mesh)
	$(CC) test.mesh.simplicialmesh.cpp $(Objects) -o test.mesh.simplicialmesh.out

test.mesh.manifold.2D.out: test.mesh.manifold.2D.cpp $(Includes.mesh)
	$(CC) test.mesh.manifold.2D.cpp $(Objects) -o test.mesh.manifold.2D.out

test.mesh.vtkwriter.out: test.mesh.vtkwriter.cpp $(Includes.mesh)
	$(CC) test.mesh.vtkwriter.cpp $(Objects) -o test.mesh.vtkwriter.out

test.mesh.vtkwriter.manifold2D.out: test.mesh.vtkwriter.manifold2D.cpp $(Includes.mesh)
	$(CC) test.mesh.vtkwriter.manifold2D.cpp $(Objects) -o test.mesh.vtkwriter.manifold2D.out

test.mesh.io.out: test.mesh.io.cpp $(Includes.mesh)
	$(CC) test.mesh.io.cpp $(Objects) -o test.mesh.io.out

mesh.tests = test.mesh.coordinates.out test.mesh.simplicialmesh.out test.mesh.manifold.2D.out test.mesh.io.out test.mesh.vtkwriter.out test.mesh.vtkwriter.manifold2D.out



test.container.out: test.container.cpp 
	$(CC) test.container.cpp  -o test.container.out


    
all: $(combinatorics.tests) $(operators.tests) $(dense.tests) $(mesh.tests) 

# $(sparse.tests) $(solver.tests)

include ../makefile.clean


