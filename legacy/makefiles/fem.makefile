# TODO: This makefile still needs to be written properly

SHELL = /bin/sh

default: all

include ../common.make 



# global.lagrangemassmatrix.o: global.lagrangemassmatrix.hpp global.lagrangemassmatrix.cpp 
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) global.lagrangemassmatrix.cpp -c -o global.lagrangemassmatrix.o 

objects = 

    
all: $(objects)  

include ../makefile.clean
