
# CC = g++ -O2 -std=c++11 -Wall -Wextra -pedantic
include ../makefile.rules 

default: all

iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
	$(CC) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.hpp crm.cpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CC) crm.cpp -c -o crm.o 

necrm.o: necrm.hpp necrm.cpp ../operators/floatvector.hpp crm.hpp 
	$(CC) necrm.cpp -c -o necrm.o 

objects = iterativesolver.o crm.o necrm.o 

unittest.crm.out: unittest.crm.cpp crm.o $(objects)
	$(CC) unittest.crm.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o unittest.crm.out

moretesting.out: moretesting.cpp crm.o $(objects)
	$(CC) moretesting.cpp ../combinatorics/*.o ../operators/*.o ../dense/*.o ../sparse/*.o *.o -o moretesting.out

unittests = unittest.crm.out moretesting.out 

all: $(objects) $(unittests)

clean:
	rm -f *.exe *.o *.gch *.exe.stackdump *.out 