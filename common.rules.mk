SHELL = /bin/bash

# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 


default: all

dirname := $(notdir $(shell pwd))
depdir  := .deps

sources      := $(wildcard *.cpp)
objects      := $(patsubst %.cpp,%.o,$(sources))
dependencies := $(patsubst %.cpp,.deps/%.d,$(sources))

sharedlibrarybasename := lib$(dirname)
sharedlibrary         := lib$(dirname).so


$(depdir): ; @mkdir -p $@

.PHONY: make_dependencies
make_dependencies: $(depdir)
	for item in $(sources); do g++ -MM $$item -MF .deps/$*.d; done

$(objects): %.o: %.cpp | $(depdir)
	@g++ -MM $*.cpp -MF .deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c -o $@ 

# NOTE: Original recipe for the shared library, now there is just one object
# $(sharedlibrary): $(objects)
# 	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

$(sharedlibrarybasename).o: $(sources) .all.cpp | $(depdir)
	@g++ -MM .all.cpp -MF .deps/.all.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) .all.cpp -c -o $@ 

$(sharedlibrary): $(sharedlibrarybasename).o
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)






-include $(dependencies)

*.o: ../common.recipe.mk ../common.rules.mk ../common.upkeep.mk ./makefile


buildobjects: $(objects)
buildso: $(sharedlibrary)

all: buildso
#buildobjects # NOTE: the .o files were required for our .so files originally


.PHONY:  list_of_objects
.SILENT: list_of_objects
list_of_objects: 
	@echo $(dirname);
	@echo $(sharedlibrary);
	@echo $(objects);
	@echo $(dependencies);


# $(sharedlibrary): #$(objects)
# 	# TODO : create the dependencies of the cpp file 
# 	MYFILE=''; for mycpp in *.cpp; do MYFILE+=$$'#include \"'$$mycpp$$'\"\n'; done; (echo "$$MYFILE") > $(sharedlibrarybasename).code
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -x c++ $(sharedlibrarybasename).code -c -o $(sharedlibrarybasename).o
# 	$(CXX) $(CXXFLAGS) -shared -o $@ $(sharedlibrarybasename).o $(LDLIBS)

