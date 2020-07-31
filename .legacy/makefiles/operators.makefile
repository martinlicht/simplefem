SHELL = /bin/sh

default: all

objects = floatvector.o io.floatvector.o linearoperator.o productoperator.o sumoperator.o scalingoperator.o diagonaloperator.o blockdiagonaloperator.o blockoperator.o 

all: $(objects) liboperators.so

include ../common.make 

include ../makefile.clean

liboperators.so: $(objects)
	g++ -shared -o liboperators.so $(objects)

floatvector.o: floatvector.cpp floatvector.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) floatvector.cpp -c -o floatvector.o 

io.floatvector.o: io.floatvector.cpp io.floatvector.hpp floatvector.hpp floatvector.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) io.floatvector.cpp -c -o io.floatvector.o 

linearoperator.o: linearoperator.cpp linearoperator.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) linearoperator.cpp -c -o linearoperator.o 

productoperator.o: productoperator.cpp productoperator.hpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) productoperator.cpp -c -o productoperator.o 

sumoperator.o: sumoperator.cpp sumoperator.hpp floatvector.hpp linearoperator.hpp scalingoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) sumoperator.cpp -c -o sumoperator.o 

scalingoperator.o: scalingoperator.cpp scalingoperator.hpp floatvector.hpp linearoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) scalingoperator.cpp -c -o scalingoperator.o 

diagonaloperator.o: diagonaloperator.cpp diagonaloperator.hpp floatvector.hpp linearoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) diagonaloperator.cpp -c -o diagonaloperator.o 

blockdiagonaloperator.o: blockdiagonaloperator.cpp blockdiagonaloperator.hpp floatvector.hpp linearoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) blockdiagonaloperator.cpp -c -o blockdiagonaloperator.o 

blockoperator.o: blockoperator.cpp blockoperator.hpp floatvector.hpp linearoperator.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) blockoperator.cpp -c -o blockoperator.o 
