SHELL = /bin/sh

default: all

objects = iterativesolver.o crm.o necrm.o 

all: $(objects) libsolver.so

include ../makefile.rules 

include ../makefile.clean

libsolver.so: $(objects)
	g++ -shared -o libsolver.so $(objects)

iterativesolver.o: iterativesolver.hpp iterativesolver.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.hpp crm.cpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) crm.cpp -c -o crm.o 

descent.residual.o: descent.residual.hpp descent.residual.cpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) descent.residual.cpp -c -o descent.residual.o 

necrm.o: necrm.hpp necrm.cpp ../operators/floatvector.hpp crm.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) necrm.cpp -c -o necrm.o 
