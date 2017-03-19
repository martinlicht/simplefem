
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

floatvector.o: floatvector.hpp floatvector.cpp 
	$(CC) floatvector.cpp -c -o floatvector.o 

linearoperator.o: linearoperator.hpp linearoperator.cpp 
	$(CC) linearoperator.cpp -c -o linearoperator.o 

productoperator.o: productoperator.hpp productoperator.cpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CC) productoperator.cpp -c -o productoperator.o 

sumoperator.o: sumoperator.hpp sumoperator.cpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CC) sumoperator.cpp -c -o sumoperator.o 

scalingoperator.o: scalingoperator.hpp scalingoperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) scalingoperator.cpp -c -o scalingoperator.o 

diagonaloperator.o: diagonaloperator.hpp diagonaloperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) diagonaloperator.cpp -c -o diagonaloperator.o 

blockdiagonaloperator.o: blockdiagonaloperator.hpp blockdiagonaloperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) blockdiagonaloperator.cpp -c -o blockdiagonaloperator.o 

# # densematrix.o: densematrix.hpp densematrix.cpp floatvector.hpp linearoperator.hpp 
# # 	$(CC) densematrix.cpp -c -o densematrix.o 
# 
# sparsematrix.o: sparsematrix.hpp sparsematrix.cpp floatvector.hpp linearoperator.hpp 
# 	$(CC) sparsematrix.cpp -c -o sparsematrix.o 
# 
# # iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
# # 	$(CC) iterativesolver.cpp -c -o iterativesolver.o 
# # 
# # crm.o: crm.hpp crm.cpp floatvector.hpp iterativesolver.hpp 
# # 	$(CC) crm.cpp -c -o crm.o 
# # 
# # necrm.o: necrm.hpp necrm.cpp floatvector.hpp densematrix.hpp crm.hpp 
# # 	$(CC) necrm.cpp -c -o necrm.o 

# matrixalgorithm.o: matrixalgorithm.cpp densematrix.hpp floatvector.hpp crm.hpp 
# 	$(CC) matrixalgorithm.cpp -c -o matrixalgorithm.o 

objects = floatvector.o linearoperator.o productoperator.o sumoperator.o scalingoperator.o diagonaloperator.o blockdiagonaloperator.o 

# densematrix.o matrixalgorithm.o sparsematrix.o iterativesolver.o crm.o necrm.o 


unittest.floatvector.out: floatvector.o floatvector.hpp unittest.floatvector.cpp $(objects) 
	$(CC) unittest.floatvector.cpp floatvector.o -o unittest.floatvector.out

unittest.scalingoperator.out: unittest.scalingoperator.cpp scalingoperator.o $(objects)
	$(CC) unittest.scalingoperator.cpp scalingoperator.o floatvector.o linearoperator.o -o unittest.scalingoperator.out
	
unittest.diagonaloperator.out: unittest.diagonaloperator.cpp diagonaloperator.o $(objects)
	$(CC) unittest.diagonaloperator.cpp diagonaloperator.o scalingoperator.o floatvector.o linearoperator.o -o unittest.diagonaloperator.out
	
# unittest.sparsematrix.out: unittest.sparsematrix.cpp sparsematrix.o $(objects)
# 	$(CC) unittest.sparsematrix.cpp ../combinatorics/*.o *.o -o unittest.sparsematrix.out
# 
# unittest.crm.out: unittest.crm.cpp crm.o $(objects)
# 	$(CC) unittest.crm.cpp ../combinatorics/*.o *.o -o unittest.crm.out


# # unittest.densematrix.out: unittest.densematrix.cpp densematrix.o $(objects) 
# # 	$(CC) unittest.densematrix.cpp ../combinatorics/*.o *.o -o unittest.densematrix.out 
# # 
# # unittest.matrixalgorithm.out: unittest.matrixalgorithm.cpp matrixalgorithm.o $(objects) 
# # 	$(CC) unittest.matrixalgorithm.cpp ../combinatorics/*.o *.o -o unittest.matrixalgorithm.out

unittests = unittest.floatvector.out unittest.scalingoperator.out unittest.diagonaloperator.out 

# # unittest.densematrix.out unittest.matrixalgorithm.out unittest.sparsematrix.out unittest.crm.out moretesting.out
    
all: $(objects) $(unittests)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 