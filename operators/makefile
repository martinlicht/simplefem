
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

floatvector.o: floatvector.hpp floatvector.cpp 
	$(CC) floatvector.cpp -c -o floatvector.o 

unittest.floatvector: floatvector.o floatvector.hpp unittest.floatvector.cpp 
	$(CC) unittest.floatvector.cpp floatvector.o -o unittest.floatvector

linearoperator.o: linearoperator.hpp linearoperator.cpp 
	$(CC) linearoperator.cpp -c -o linearoperator.o 

scalingoperator.o: scalingoperator.hpp scalingoperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) scalingoperator.cpp -c -o scalingoperator.o 

unittest.scalingoperator: scalingoperator.o floatvector.o unittest.scalingoperator.cpp 
	$(CC) unittest.scalingoperator.cpp scalingoperator.o floatvector.o linearoperator.o -o unittest.scalingoperator
	
sparsematrix.o: sparsematrix.hpp sparsematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) sparsematrix.cpp -c -o sparsematrix.o 

unittest.sparsematrix: sparsematrix.o floatvector.o unittest.sparsematrix.cpp 
	$(CC) unittest.sparsematrix.cpp sparsematrix.o scalingoperator.o floatvector.o linearoperator.o -o unittest.sparsematrix

iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
	$(CC) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.hpp crm.cpp floatvector.hpp iterativesolver.hpp iterativesolver.o 
	$(CC) crm.cpp -c -o crm.o 

unittest.crm: crm.o floatvector.o unittest.crm.cpp 
	$(CC) unittest.crm.cpp crm.o iterativesolver.o scalingoperator.o floatvector.o linearoperator.o sparsematrix.o -o unittest.crm

# moretesting: crm.o floatvector.o moretesting.cpp 
# 	$(CC) moretesting.cpp crm.o iterativesolver.o scalingoperator.o floatvector.o linearoperator.o sparsematrix.o -o moretesting
moretesting: crm.o floatvector.o moretesting.cpp 
	$(CC) moretesting.cpp -o moretesting

densematrix.o: densematrix.hpp densematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

unittest.densematrix: densematrix.o floatvector.o linearoperator.o unittest.densematrix.cpp 
	$(CC) unittest.densematrix.cpp densematrix.o ../combinatorics/natmap.o floatvector.o linearoperator.o -o unittest.densematrix 

matrixalgorithm.o: matrixalgorithm.cpp densematrix.hpp floatvector.hpp crm.hpp 
	$(CC) matrixalgorithm.cpp -c -o matrixalgorithm.o 

unittest.matrixalgorithm: matrixalgorithm.o unittest.matrixalgorithm.cpp 
	$(CC) unittest.matrixalgorithm.cpp densematrix.o matrixalgorithm.o crm.o iterativesolver.o scalingoperator.o ../combinatorics/natmap.o floatvector.o linearoperator.o -o unittest.matrixalgorithm

all: floatvector.o unittest.floatvector linearoperator.o scalingoperator.o unittest.scalingoperator sparsematrix.o unittest.sparsematrix crm.o unittest.crm densematrix.o unittest.densematrix moretesting matrixalgorithm.o unittest.matrixalgorithm

clean: 
	rm *.exe *.o *.gch *.exe.stackdump 