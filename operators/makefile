
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

floatvector.o: floatvector.hpp floatvector.cpp 
	$(CC) floatvector.cpp -c -o floatvector.o 

unittest.floatvector.out: floatvector.o floatvector.hpp unittest.floatvector.cpp 
	$(CC) unittest.floatvector.cpp floatvector.o -o unittest.floatvector.out

linearoperator.o: linearoperator.hpp linearoperator.cpp 
	$(CC) linearoperator.cpp -c -o linearoperator.o 

scalingoperator.o: scalingoperator.hpp scalingoperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) scalingoperator.cpp -c -o scalingoperator.o 

productoperator.o: productoperator.hpp productoperator.cpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CC) productoperator.cpp -c -o productoperator.o 

sumoperator.o: sumoperator.hpp sumoperator.cpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CC) sumoperator.cpp -c -o sumoperator.o 

unittest.scalingoperator.out: scalingoperator.o floatvector.o unittest.scalingoperator.cpp 
	$(CC) unittest.scalingoperator.cpp scalingoperator.o floatvector.o linearoperator.o -o unittest.scalingoperator.out
	
diagonaloperator.o: diagonaloperator.hpp diagonaloperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) diagonaloperator.cpp -c -o diagonaloperator.o 

unittest.diagonaloperator.out: diagonaloperator.o floatvector.o unittest.diagonaloperator.cpp 
	$(CC) unittest.diagonaloperator.cpp diagonaloperator.o floatvector.o linearoperator.o -o unittest.diagonaloperator.out
	
sparsematrix.o: sparsematrix.hpp sparsematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) sparsematrix.cpp -c -o sparsematrix.o 

unittest.sparsematrix.out: sparsematrix.o floatvector.o unittest.sparsematrix.cpp 
	$(CC) unittest.sparsematrix.cpp sparsematrix.o scalingoperator.o floatvector.o linearoperator.o -o unittest.sparsematrix.out

iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
	$(CC) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.hpp crm.cpp floatvector.hpp iterativesolver.hpp iterativesolver.o 
	$(CC) crm.cpp -c -o crm.o 

unittest.crm.out: crm.o floatvector.o unittest.crm.cpp 
	$(CC) unittest.crm.cpp crm.o iterativesolver.o scalingoperator.o floatvector.o linearoperator.o sparsematrix.o -o unittest.crm.out

# moretesting: crm.o floatvector.o moretesting.cpp 
# 	$(CC) moretesting.cpp crm.o iterativesolver.o scalingoperator.o floatvector.o linearoperator.o sparsematrix.o -o moretesting
moretesting: crm.o floatvector.o moretesting.cpp 
	$(CC) moretesting.cpp -o moretesting

densematrix.o: densematrix.hpp densematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

unittest.densematrix.out: densematrix.o floatvector.o linearoperator.o unittest.densematrix.cpp 
	$(CC) unittest.densematrix.cpp densematrix.o ../combinatorics/indexrange.o ../combinatorics/indexmap.o floatvector.o linearoperator.o matrixalgorithm.o productoperator.o -o unittest.densematrix.out 

matrixalgorithm.o: matrixalgorithm.cpp densematrix.hpp floatvector.hpp crm.hpp 
	$(CC) matrixalgorithm.cpp -c -o matrixalgorithm.o 

unittest.matrixalgorithm.out: matrixalgorithm.o unittest.matrixalgorithm.cpp 
	$(CC) unittest.matrixalgorithm.cpp densematrix.o matrixalgorithm.o crm.o iterativesolver.o scalingoperator.o ../combinatorics/indexrange.o ../combinatorics/indexmap.o floatvector.o linearoperator.o -o unittest.matrixalgorithm.out

all: floatvector.o unittest.floatvector.out linearoperator.o scalingoperator.o unittest.scalingoperator.out diagonaloperator.o unittest.diagonaloperator.out sparsematrix.o unittest.sparsematrix.out crm.o unittest.crm.out densematrix.o unittest.densematrix.out moretesting matrixalgorithm.o unittest.matrixalgorithm.out sumoperator.o productoperator.o

clean:
	rm *.exe *.o *.gch *.exe.stackdump *.out 