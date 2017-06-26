
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

floatvector.o: floatvector.hpp floatvector.cpp 
	$(CC) floatvector.cpp -c -o floatvector.o 

io.floatvector.o: io.floatvector.hpp io.floatvector.cpp floatvector.hpp floatvector.cpp 
	$(CC) io.floatvector.cpp -c -o io.floatvector.o 

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

blockoperator.o: blockoperator.hpp blockoperator.cpp floatvector.hpp linearoperator.hpp
	$(CC) blockoperator.cpp -c -o blockoperator.o 



objects = floatvector.o io.floatvector.o linearoperator.o productoperator.o sumoperator.o scalingoperator.o diagonaloperator.o blockdiagonaloperator.o blockoperator.o 



all: $(objects) 

include ../makefile.clean

