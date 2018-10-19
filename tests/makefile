SHELL = /bin/sh

include ../makefile.rules 

default: all

Objects.basic = 
Objects.combinatorics = ../combinatorics/*.o
Objects.operators = ../operators/*.o
Objects.dense = ../dense/*.o
Objects.sparse = ../sparse/*.o
Objects.solver = ../solver/*.o
Objects.mesh = ../mesh/*.o
Objects.vtk = ../vtk/*.o
Objects.matrixmarket = ../matrixmarket/*.o

Headers.basic = ../basic/*.hpp
Headers.combinatorics = ../combinatorics/*.hpp
Headers.operators = ../operators/*.hpp
Headers.dense = ../dense/*.hpp
Headers.sparse = ../sparse/*.hpp
Headers.solver = ../solver/*.hpp
Headers.mesh = ../mesh/*.hpp
Headers.vtk = ../vtk/*.hpp
Headers.matrixmarket = ../matrixmarket/*.hpp

Includes.basic = $(Headers.basic) $(Objects.basic) 
Includes.combinatorics = $(Headers.combinatorics) $(Objects.combinatorics) 
Includes.operators = $(Headers.operators) $(Objects.operators) $(Includes.combinatorics)
Includes.dense = $(Headers.dense) $(Objects.dense) $(Includes.operators)
Includes.sparse = $(Headers.sparse) $(Objects.sparse) $(Includes.operators)
Includes.solver = $(Headers.solver) $(Objects.solver) $(Includes.operators)
Includes.mesh = $(Headers.mesh) $(Objects.mesh) $(Objects.combinatorics) $(Includes.operators)
Includes.vtk = $(Headers.vtk) $(Objects.vtk) $(Includes.mesh)
Includes.matrixmarket = $(Headers.matrixmarket) $(Objects.matrixmarket)

Links.basic         = $(Objects.basic) 
Links.combinatorics = $(Objects.combinatorics) 
Links.operators     = $(Objects.operators) $(Objects.combinatorics)
Links.dense         = $(Objects.dense) $(Objects.sparse) $(Objects.operators) $(Objects.combinatorics)
Links.sparse        = $(Objects.sparse) $(Objects.operators)
Links.solver        = $(Objects.solver) $(Objects.operators)
Links.mesh          = $(Objects.mesh) $(Objects.combinatorics) $(Objects.operators) $(Objects.dense) $(Objects.sparse) 
Links.vtk           = $(Objects.vtk) $(Objects.mesh) $(Objects.dense) $(Objects.sparse) $(Objects.combinatorics) $(Objects.operators)
Links.matrixmarket  = $(Objects.matrixmarket) 


Objects = ../basic/*.o ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o ../solver/*.o ../mesh/*.o ../vtk/*.o ../matrixmarket/*.o
Headers = ../basic/*.hpp ../combinatorics/*.hpp ../operators/*.hpp ../dense/*.hpp ../sparse/*.hpp ../solver/*.hpp ../mesh/*.hpp ../vtk/*.hpp ../matrixmarket/*.hpp

Includes = $(Headers) $(Objects) 




basic/sorthack.out: basic/sorthack.cpp $(Includes.basic)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) basic/sorthack.cpp $(Links.basic) -o basic/sorthack.out 

basic/logging.out: basic/logging.cpp $(Includes.basic)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) basic/logging.cpp $(Links.basic) -o basic/logging.out 

basic.tests = basic/sorthack.out basic/logging.out





combinatorics/indexrange.out: combinatorics/indexrange.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/indexrange.cpp $(Links.combinatorics) -o combinatorics/indexrange.out

combinatorics/indexmap.out: combinatorics/indexmap.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/indexmap.cpp $(Links.combinatorics) -o combinatorics/indexmap.out

combinatorics/multiindex.out: combinatorics/multiindex.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/multiindex.cpp $(Links.combinatorics) -o combinatorics/multiindex.out

combinatorics/generateindexmaps.out: combinatorics/generateindexmaps.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/generateindexmaps.cpp $(Links.combinatorics) -o combinatorics/generateindexmaps.out

combinatorics/generatemultiindices.out: combinatorics/generatemultiindices.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/generatemultiindices.cpp $(Links.combinatorics) -o combinatorics/generatemultiindices.out

combinatorics/heappermgen.out: combinatorics/heappermgen.cpp $(Includes.combinatorics) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) combinatorics/heappermgen.cpp $(Links.combinatorics) -o combinatorics/heappermgen.out

combinatorics.tests = combinatorics/indexrange.out combinatorics/indexmap.out combinatorics/multiindex.out combinatorics/generateindexmaps.out combinatorics/generatemultiindices.out combinatorics/heappermgen.out 
    



operators/floatvector.out: operators/floatvector.cpp $(Includes.operators) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/floatvector.cpp $(Links.operators) -o operators/floatvector.out

operators/scalingoperator.out: operators/scalingoperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/scalingoperator.cpp $(Links.operators) -o operators/scalingoperator.out
	
operators/diagonaloperator.out: operators/diagonaloperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/diagonaloperator.cpp $(Links.operators) -o operators/diagonaloperator.out
	
operators/sumoperator.out: operators/sumoperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/sumoperator.cpp $(Links.operators) -o operators/sumoperator.out
	
operators/productoperator.out: operators/productoperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/productoperator.cpp $(Links.operators) -o operators/productoperator.out
	
operators/blockdiagonaloperator.out: operators/blockdiagonaloperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/blockdiagonaloperator.cpp $(Links.operators) -o operators/blockdiagonaloperator.out
	
operators/blockoperator.out: operators/blockoperator.cpp $(Includes.operators)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) operators/blockoperator.cpp $(Links.operators) -o operators/blockoperator.out
	
operators.tests = operators/floatvector.out operators/scalingoperator.out operators/diagonaloperator.out operators/sumoperator.out operators/productoperator.out operators/blockdiagonaloperator.out operators/blockoperator.out


    
    
    
    
dense/densematrix.out: dense/densematrix.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/densematrix.cpp $(Links.dense) -o dense/densematrix.out 

dense/transpose.out: dense/transpose.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/transpose.cpp $(Links.dense) -o dense/transpose.out

dense/scalarfunctions.out: dense/scalarfunctions.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/scalarfunctions.cpp $(Links.dense) -o dense/scalarfunctions.out

dense/detinvcof.out: dense/detinvcof.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/detinvcof.cpp $(Links.dense) -o dense/detinvcof.out

dense/matrixalgorithm.out: dense/matrixalgorithm.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/matrixalgorithm.cpp $(Links.dense) -o dense/matrixalgorithm.out

dense/cholesky.out: dense/cholesky.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/cholesky.cpp $(Links.dense) -o dense/cholesky.out

dense/iterativeinverse.out: dense/iterativeinverse.cpp $(Includes.dense) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) dense/iterativeinverse.cpp $(Links.dense) -o dense/iterativeinverse.out

dense.tests = dense/densematrix.out dense/transpose.out dense/detinvcof.out dense/matrixalgorithm.out dense/cholesky.out dense/scalarfunctions.out dense/iterativeinverse.out




sparse/sparsematrix.out: sparse/sparsematrix.cpp $(Includes.sparse)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) sparse/sparsematrix.cpp $(Links.sparse) -o sparse/sparsematrix.out

sparse.tests = sparse/sparsematrix.out 






solver/crm.out: solver/crm.cpp $(Includes.solver)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) solver/crm.cpp $(Links.solver) -o solver/crm.out

solver/moretesting.out: solver/moretesting.cpp $(Includes.solver)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) solver/moretesting.cpp $(Links.solver) -o solver/moretesting.out

solver.tests = solver/crm.out solver/moretesting.out 





mesh/coordinates.out: mesh/coordinates.cpp $(Includes.mesh) 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/coordinates.cpp $(Links.mesh) -o mesh/coordinates.out

# mesh/simplicialmesh.out: mesh/simplicialmesh.cpp $(Includes.mesh)
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/simplicialmesh.cpp $(Links.mesh) -o mesh/simplicialmesh.out

mesh/simplicial1D.out: mesh/simplicial1D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/simplicial1D.cpp $(Links.mesh) -o mesh/simplicial1D.out

mesh/io.simplicial1D.out: mesh/io.simplicial1D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/io.simplicial1D.cpp $(Links.mesh) -o mesh/io.simplicial1D.out

mesh/uniformrefinement1D.out: mesh/uniformrefinement1D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/uniformrefinement1D.cpp $(Links.mesh) -o mesh/uniformrefinement1D.out

# mesh/manifold2D.out: mesh/manifold2D.cpp $(Includes.mesh)
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/manifold2D.cpp $(Links.mesh) -o mesh/manifold2D.out

mesh/simplicial2D.out: mesh/simplicial2D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/simplicial2D.cpp $(Links.mesh) -o mesh/simplicial2D.out

mesh/io.simplicial2D.out: mesh/io.simplicial2D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/io.simplicial2D.cpp $(Links.mesh) -o mesh/io.simplicial2D.out

mesh/2DBisection.out: mesh/2DBisection.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/2DBisection.cpp $(Links.mesh) -o mesh/2DBisection.out

mesh/2Dleb.out: mesh/2Dleb.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/2Dleb.cpp $(Links.mesh) -o mesh/2Dleb.out

mesh/2Dnvb.out: mesh/2Dnvb.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/2Dnvb.cpp $(Links.mesh) -o mesh/2Dnvb.out

mesh/midpointrefine.out: mesh/midpointrefine.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/midpointrefine.cpp $(Links.mesh) -o mesh/midpointrefine.out

mesh/simplicial3D.out: mesh/simplicial3D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/simplicial3D.cpp $(Links.mesh) -o mesh/simplicial3D.out

mesh/io.simplicial3D.out: mesh/io.simplicial3D.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/io.simplicial3D.cpp $(Links.mesh) -o mesh/io.simplicial3D.out

mesh/3DBisection.out: mesh/3DBisection.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/3DBisection.cpp $(Links.mesh) -o mesh/3DBisection.out

mesh/3Dleb.out: mesh/3Dleb.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/3Dleb.cpp $(Links.mesh) -o mesh/3Dleb.out

mesh/simplicialND.out: mesh/simplicialND.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/simplicialND.cpp $(Links.mesh) -o mesh/simplicialND.out

mesh/io.simplicialND.out: mesh/io.simplicialND.cpp $(Includes.mesh)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/io.simplicialND.cpp $(Links.mesh) -o mesh/io.simplicialND.out

# mesh/mesh.manifold.2D.out: mesh/mesh.manifold.2D.cpp $(Includes.mesh)
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/mesh.manifold.2D.cpp $(Links.mesh) -o mesh/mesh.manifold.2D.out

# mesh/io.out: mesh/io.cpp $(Includes.mesh)
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) mesh/io.cpp $(Links.mesh) -o mesh/io.out

mesh.tests = mesh/coordinates.out mesh/simplicial1D.out mesh/io.simplicial1D.out mesh/uniformrefinement1D.out mesh/simplicial2D.out mesh/2DBisection.out mesh/2Dleb.out mesh/2Dnvb.out mesh/midpointrefine.out mesh/io.simplicial2D.out mesh/simplicial3D.out mesh/3DBisection.out mesh/3Dleb.out mesh/io.simplicial3D.out mesh/simplicialND.out mesh/io.simplicialND.out
# mesh/simplicial3D.out mesh/io.out mesh/vtkwriter.out mesh/mesh.manifold.2D.out  mesh/manifold2D.out





# vtk/vtkwriter.out: vtk/vtkwriter.cpp $(Includes.vtk)
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.cpp $(Links.vtk) -o vtk/vtkwriter.out

vtk/vtkwriter.mesh1D.out: vtk/vtkwriter.mesh1D.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh1D.cpp $(Links.vtk) -o vtk/vtkwriter.mesh1D.out

vtk/vtkwriter.mesh2D.out: vtk/vtkwriter.mesh2D.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh2D.cpp $(Links.vtk) -o vtk/vtkwriter.mesh2D.out

vtk/vtkwriter.mesh3D.out: vtk/vtkwriter.mesh3D.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh3D.cpp $(Links.vtk) -o vtk/vtkwriter.mesh3D.out

vtk/vtkwriter.mesh2Dnvb.out: vtk/vtkwriter.mesh2Dnvb.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh2Dnvb.cpp $(Links.vtk) -o vtk/vtkwriter.mesh2Dnvb.out

vtk/vtkwriter.mesh2Dleb.out: vtk/vtkwriter.mesh2Dleb.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh2Dleb.cpp $(Links.vtk) -o vtk/vtkwriter.mesh2Dleb.out

vtk/vtkwriter.mesh2Dur.out: vtk/vtkwriter.mesh2Dur.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh2Dur.cpp $(Links.vtk) -o vtk/vtkwriter.mesh2Dur.out

vtk/vtkwriter.mesh3Dur.out: vtk/vtkwriter.mesh3Dur.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh3Dur.cpp $(Links.vtk) -o vtk/vtkwriter.mesh3Dur.out

vtk/vtkwriter.mesh3Dleb.out: vtk/vtkwriter.mesh3Dleb.cpp $(Includes.vtk)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) vtk/vtkwriter.mesh3Dleb.cpp $(Links.vtk) -o vtk/vtkwriter.mesh3Dleb.out

vtk.tests = vtk/vtkwriter.mesh1D.out vtk/vtkwriter.mesh2D.out vtk/vtkwriter.mesh2Dnvb.out vtk/vtkwriter.mesh2Dleb.out vtk/vtkwriter.mesh2Dur.out vtk/vtkwriter.mesh3D.out vtk/vtkwriter.mesh3Dur.out vtk/vtkwriter.mesh3Dleb.out



matrixmarket/mm.out: matrixmarket/mm.cpp $(Includes.matrixmarket)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) matrixmarket/mm.cpp $(Links.matrixmarket) -o matrixmarket/mm.out

matrixmarket.tests = matrixmarket/mm.out






test.compiler.out: test.compiler.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test.compiler.cpp  -o test.compiler.out

test.container.out: test.container.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test.container.cpp  -o test.container.out



all: test.compiler.out $(basic.tests) $(combinatorics.tests) $(operators.tests) $(dense.tests) $(mesh.tests) $(vtk.tests) $(matrixmarket.tests) 

# $(sparse.tests) $(solver.tests)


include makefile.tests


clean: 
	-rm -f *.o *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 
	-rm -f */*.o */*.gch
	-rm -f */*.exe */*.exe.stackdump
	-rm -f */*.out */*.out.stackdump 


