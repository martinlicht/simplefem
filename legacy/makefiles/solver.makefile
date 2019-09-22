SHELL = /bin/sh

default: all

objects = iterativesolver.o crm.o pcrm.o necrm.o 

all: $(objects) libsolver.so

include ../common.make 

include ../makefile.clean

libsolver.so: $(objects)
	g++ -shared -o libsolver.so $(objects)

iterativesolver.o: iterativesolver.cpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) iterativesolver.cpp -c -o iterativesolver.o 

crm.o: crm.cpp crm.hpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) crm.cpp -c -o crm.o 

pcrm.o: pcrm.cpp pcrm.hpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) pcrm.cpp -c -o pcrm.o 

descent.residual.o: descent.residual.hpp descent.residual.cpp ../operators/floatvector.hpp iterativesolver.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) descent.residual.cpp -c -o descent.residual.o 

necrm.o: necrm.cpp necrm.hpp ../operators/floatvector.hpp crm.hpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) necrm.cpp -c -o necrm.o 
