
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

floatvector.o: floatvector.hpp floatvector.cpp 
	$(CC) floatvector.cpp -c -o floatvector.o 

unittest_floatvector: floatvector.o floatvector.hpp unittest_floatvector.cpp 
	$(CC) unittest_floatvector.cpp floatvector.o -o unittest_floatvector

linearoperator.o: linearoperator.hpp linearoperator.cpp 
	$(CC) linearoperator.cpp -c -o linearoperator.o 

scalingoperator.o: scalingoperator.hpp scalingoperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) scalingoperator.cpp -c -o scalingoperator.o 

unittest_scalingoperator: scalingoperator.o floatvector.o unittest_scalingoperator.cpp 
	$(CC) unittest_scalingoperator.cpp scalingoperator.o floatvector.o linearoperator.o -o unittest_scalingoperator
	
sparsematrix.o: sparsematrix.hpp sparsematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) sparsematrix.cpp -c -o sparsematrix.o 

unittest_sparsematrix: sparsematrix.o floatvector.o unittest_sparsematrix.cpp 
	$(CC) unittest_sparsematrix.cpp sparsematrix.o scalingoperator.o floatvector.o linearoperator.o -o unittest_sparsematrix

iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
	$(CC) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.hpp crm.cpp floatvector.hpp iterativesolver.hpp iterativesolver.o 
	$(CC) crm.cpp -c -o crm.o 

unittest_crm: crm.o floatvector.o unittest_crm.cpp 
	$(CC) unittest_crm.cpp crm.o iterativesolver.o scalingoperator.o floatvector.o linearoperator.o sparsematrix.o -o unittest_crm

densematrix.o: densematrix.hpp densematrix.cpp floatvector.hpp linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

unittest_densematrix: densematrix.o floatvector.o linearoperator.o unittest_densematrix.cpp 
	$(CC) unittest_densematrix.cpp densematrix.o ../combinatorics/natmap.o floatvector.o linearoperator.o -o unittest_densematrix 

matrixalgorithm.o: matrixalgorithm.cpp densematrix.hpp floatvector.hpp crm.hpp 
	$(CC) matrixalgorithm.cpp -c -o matrixalgorithm.o 

unittest_matrixalgorithm: matrixalgorithm.o unittest_matrixalgorithm.cpp 
	$(CC) unittest_matrixalgorithm.cpp densematrix.o matrixalgorithm.o crm.o iterativesolver.o scalingoperator.o ../combinatorics/natmap.o floatvector.o linearoperator.o -o unittest_matrixalgorithm

all: floatvector.o unittest_floatvector linearoperator.o scalingoperator.o unittest_scalingoperator sparsematrix.o unittest_sparsematrix crm.o unittest_crm densematrix.o unittest_densematrix matrixalgorithm.o unittest_matrixalgorithm
