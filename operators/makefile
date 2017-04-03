
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



objects = floatvector.o linearoperator.o productoperator.o sumoperator.o scalingoperator.o diagonaloperator.o blockdiagonaloperator.o 



all: $(objects) 

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out *.out.stackdump 
