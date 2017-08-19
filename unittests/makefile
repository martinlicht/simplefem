
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

Objects.combinatorics = ../combinatorics/*.o
Objects.operators = ../operators/*.o
Objects.dense = ../dense/*.o
Objects.sparse = ../sparse/*.o
Objects.solver = ../solver/*.o
Objects.mesh = ../mesh/*.o
Objects.vtk = ../vtk/*.o

Headers.combinatorics = ../combinatorics/*.o
Headers.operators = ../operators/*.o
Headers.dense = ../dense/*.o
Headers.sparse = ../sparse/*.o
Headers.solver = ../solver/*.o
Headers.mesh = ../mesh/*.o
Headers.vtk = ../vtk/*.o

Includes.combinatorics = $(Headers.combinatorics) $(Objects.combinatorics) 
Includes.operators = $(Headers.operators) $(Objects.operators) $(Includes.combinatorics)
Includes.dense = $(Headers.dense) $(Objects.dense) $(Includes.operators)
Includes.sparse = $(Headers.sparse) $(Objects.sparse) $(Includes.operators)
Includes.solver = $(Headers.solver) $(Objects.solver) $(Includes.operators)
Includes.mesh = $(Headers.mesh) $(Objects.mesh) $(Objects.combinatorics) $(Includes.operators)
Includes.vtk = $(Headers.vtk) $(Objects.vtk) $(Includes.mesh)

Links.combinatorics = $(Objects.combinatorics) 
Links.operators     = $(Objects.operators) $(Objects.combinatorics)
Links.dense         = $(Objects.dense) $(Objects.sparse) $(Objects.operators) $(Objects.combinatorics)
Links.sparse        = $(Objects.sparse) $(Objects.operators)
Links.solver        = $(Objects.solver) $(Objects.operators)
Links.mesh          = $(Objects.mesh) $(Objects.combinatorics) $(Objects.operators) $(Objects.dense) $(Objects.sparse) 
Links.vtk           = $(Objects.vtk) $(Objects.mesh)


Objects = ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o ../solver/*.o ../mesh/*.o ../vtk/*.o
Headers = ../combinatorics/*.hpp ../operators/*.hpp ../dense/*.hpp ../sparse/*.hpp ../solver/*.hpp ../mesh/*.hpp ../vtk/*.hpp

Includes = $(Headers) $(Objects) 

test.combinatorics.indexrange.out: test.combinatorics.indexrange.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.indexrange.cpp $(Links.combinatorics) -o test.combinatorics.indexrange.out

test.combinatorics.indexmap.out: test.combinatorics.indexmap.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.indexmap.cpp $(Links.combinatorics) -o test.combinatorics.indexmap.out

test.combinatorics.multiindex.out: test.combinatorics.multiindex.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.multiindex.cpp $(Links.combinatorics) -o test.combinatorics.multiindex.out

test.combinatorics.generateindexmaps.out: test.combinatorics.generateindexmaps.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.generateindexmaps.cpp $(Links.combinatorics) -o test.combinatorics.generateindexmaps.out

test.combinatorics.generatemultiindices.out: test.combinatorics.generatemultiindices.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.generatemultiindices.cpp $(Links.combinatorics) -o test.combinatorics.generatemultiindices.out

test.combinatorics.heapsalgorithm.out: test.combinatorics.heapsalgorithm.cpp $(Includes.combinatorics) 
	$(CC) test.combinatorics.heapsalgorithm.cpp $(Links.combinatorics) -o test.combinatorics.heapsalgorithm.out

combinatorics.tests = test.combinatorics.indexrange.out test.combinatorics.indexmap.out test.combinatorics.multiindex.out test.combinatorics.generateindexmaps.out test.combinatorics.generatemultiindices.out test.combinatorics.heapsalgorithm.out 
    



test.operators.floatvector.out: test.operators.floatvector.cpp $(Includes.operators) 
	$(CC) test.operators.floatvector.cpp $(Links.operators) -o test.operators.floatvector.out

test.operators.scalingoperator.out: test.operators.scalingoperator.cpp $(Includes.operators)
	$(CC) test.operators.scalingoperator.cpp $(Links.operators) -o test.operators.scalingoperator.out
	
test.operators.diagonaloperator.out: test.operators.diagonaloperator.cpp $(Includes.operators)
	$(CC) test.operators.diagonaloperator.cpp $(Links.operators) -o test.operators.diagonaloperator.out
	
test.operators.sumoperator.out: test.operators.sumoperator.cpp $(Includes.operators)
	$(CC) test.operators.sumoperator.cpp $(Links.operators) -o test.operators.sumoperator.out
	
test.operators.productoperator.out: test.operators.productoperator.cpp $(Includes.operators)
	$(CC) test.operators.productoperator.cpp $(Links.operators) -o test.operators.productoperator.out
	
test.operators.blockdiagonaloperator.out: test.operators.blockdiagonaloperator.cpp $(Includes.operators)
	$(CC) test.operators.blockdiagonaloperator.cpp $(Links.operators) -o test.operators.blockdiagonaloperator.out
	
test.operators.blockoperator.out: test.operators.blockoperator.cpp $(Includes.operators)
	$(CC) test.operators.blockoperator.cpp $(Links.operators) -o test.operators.blockoperator.out
	
operators.tests = test.operators.floatvector.out test.operators.scalingoperator.out test.operators.diagonaloperator.out test.operators.sumoperator.out test.operators.productoperator.out test.operators.blockdiagonaloperator.out test.operators.blockoperator.out


    
    
    
    
test.dense.densematrix.out: test.dense.densematrix.cpp $(Includes.dense) 
	$(CC) test.dense.densematrix.cpp $(Links.dense) -o test.dense.densematrix.out 

test.dense.transpose.out: test.dense.transpose.cpp $(Includes.dense) 
	$(CC) test.dense.transpose.cpp $(Links.dense) -o test.dense.transpose.out

test.dense.scalarfunctions.out: test.dense.scalarfunctions.cpp $(Includes.dense) 
	$(CC) test.dense.scalarfunctions.cpp $(Links.dense) -o test.dense.scalarfunctions.out

test.dense.detinvcof.out: test.dense.detinvcof.cpp $(Includes.dense) 
	$(CC) test.dense.detinvcof.cpp $(Links.dense) -o test.dense.detinvcof.out

test.dense.matrixalgorithm.out: test.dense.matrixalgorithm.cpp $(Includes.dense) 
	$(CC) test.dense.matrixalgorithm.cpp $(Links.dense) -o test.dense.matrixalgorithm.out

test.dense.cholesky.out: test.dense.cholesky.cpp $(Includes.dense) 
	$(CC) test.dense.cholesky.cpp $(Links.dense) -o test.dense.cholesky.out

test.dense.iterativeinverse.out: test.dense.iterativeinverse.cpp $(Includes.dense) 
	$(CC) test.dense.iterativeinverse.cpp $(Links.dense) -o test.dense.iterativeinverse.out

dense.tests = test.dense.densematrix.out test.dense.transpose.out test.dense.detinvcof.out test.dense.matrixalgorithm.out test.dense.cholesky.out test.dense.scalarfunctions.out test.dense.iterativeinverse.out




test.sparse.sparsematrix.out: test.sparse.sparsematrix.cpp $(Includes.sparse)
	$(CC) test.sparse.sparsematrix.cpp $(Links.sparse) -o test.sparse.sparsematrix.out

sparse.tests = test.sparse.sparsematrix.out 






test.solver.crm.out: test.solver.crm.cpp $(Includes.solver)
	$(CC) test.solver.crm.cpp $(Links.solver) -o test.solver.crm.out

test.solver.moretesting.out: test.solver.moretesting.cpp $(Includes.solver)
	$(CC) test.solver.moretesting.cpp $(Links.solver) -o test.solver.moretesting.out

solver.tests = test.solver.crm.out test.solver.moretesting.out 





test.mesh.coordinates.out: test.mesh.coordinates.cpp $(Includes.mesh) 
	$(CC) test.mesh.coordinates.cpp $(Links.mesh) -o test.mesh.coordinates.out

# test.mesh.simplicialmesh.out: test.mesh.simplicialmesh.cpp $(Includes.mesh)
# 	$(CC) test.mesh.simplicialmesh.cpp $(Links.mesh) -o test.mesh.simplicialmesh.out

test.mesh.simplicial1D.out: test.mesh.simplicial1D.cpp $(Includes.mesh)
	$(CC) test.mesh.simplicial1D.cpp $(Links.mesh) -o test.mesh.simplicial1D.out

test.io.simplicial1D.out: test.io.simplicial1D.cpp $(Includes.mesh)
	$(CC) test.io.simplicial1D.cpp $(Links.mesh) -o test.io.simplicial1D.out

test.mesh.uniformrefinement1D.out: test.mesh.uniformrefinement1D.cpp $(Includes.mesh)
	$(CC) test.mesh.uniformrefinement1D.cpp $(Links.mesh) -o test.mesh.uniformrefinement1D.out

test.mesh.manifold2D.out: test.mesh.manifold2D.cpp $(Includes.mesh)
	$(CC) test.mesh.manifold2D.cpp $(Links.mesh) -o test.mesh.manifold2D.out

test.mesh.simplicial2D.out: test.mesh.simplicial2D.cpp $(Includes.mesh)
	$(CC) test.mesh.simplicial2D.cpp $(Links.mesh) -o test.mesh.simplicial2D.out

test.io.simplicial2D.out: test.io.simplicial2D.cpp $(Includes.mesh)
	$(CC) test.io.simplicial2D.cpp $(Links.mesh) -o test.io.simplicial2D.out

# test.mesh.mesh.manifold.2D.out: test.mesh.mesh.manifold.2D.cpp $(Includes.mesh)
# 	$(CC) test.mesh.mesh.manifold.2D.cpp $(Links.mesh) -o test.mesh.mesh.manifold.2D.out

# test.mesh.io.out: test.mesh.io.cpp $(Includes.mesh)
# 	$(CC) test.mesh.io.cpp $(Links.mesh) -o test.mesh.io.out

mesh.tests = test.mesh.coordinates.out test.mesh.simplicial1D.out test.io.simplicial1D.out test.mesh.uniformrefinement1D.out test.mesh.manifold2D.out test.mesh.simplicial2D.out test.io.simplicial2D.out 
# test.mesh.simplicialmesh.out test.mesh.io.out test.mesh.vtkwriter.out test.mesh.mesh.manifold.2D.out 





# test.vtk.vtkwriter.out: test.vtk.vtkwriter.cpp $(Includes.vtk)
# 	$(CC) test.vtk.vtkwriter.cpp $(Links.vtk) -o test.vtk.vtkwriter.out

test.vtk.vtkwriter.manifold2D.out: test.vtk.vtkwriter.manifold2D.cpp $(Includes.vtk)
	$(CC) test.vtk.vtkwriter.manifold2D.cpp $(Links.vtk) -o test.vtk.vtkwriter.manifold2D.out

vtk.tests = test.mesh.vtkwriter.manifold2D.out







test.container.out: test.container.cpp 
	$(CC) test.container.cpp  -o test.container.out



all: $(combinatorics.tests) $(operators.tests) $(dense.tests) $(mesh.tests) 

# $(sparse.tests) $(solver.tests)

include ../makefile.clean


