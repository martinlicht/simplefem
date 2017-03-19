
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

densematrix.o: densematrix.hpp densematrix.cpp ../operators/floatvector.hpp ../operators/linearoperator.hpp 
	$(CC) densematrix.cpp -c -o densematrix.o 

matrixalgorithm.o: matrixalgorithm.cpp densematrix.hpp ../operators/floatvector.hpp
	$(CC) matrixalgorithm.cpp -c -o matrixalgorithm.o 

objects = densematrix.o matrixalgorithm.o




# unittest.densematrix.out: unittest.densematrix.cpp densematrix.o $(objects) 
# 	$(CC) unittest.densematrix.cpp ../combinatorics/*.o ../operators/*.o *.o -o unittest.densematrix.out 
# 
# unittest.matrixalgorithm.out: unittest.matrixalgorithm.cpp matrixalgorithm.o $(objects) 
# 	$(CC) unittest.matrixalgorithm.cpp ../combinatorics/*.o ../operators/*.o *.o -o unittest.matrixalgorithm.out
# 
# unittests = unittest.densematrix.out unittest.matrixalgorithm.out
    
# all: $(objects) $(unittests)
all: $(objects)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 