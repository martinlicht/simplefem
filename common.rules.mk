SHELL = /bin/bash

# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 



dirname := $(notdir $(shell pwd))
depdir  := .deps

sources      := $(wildcard *.cpp)
headers      := $(wildcard *.hpp)
objects      := $(patsubst %.cpp,%.o,$(sources))
dependencies := $(patsubst %.cpp,.deps/%.d,$(sources))

sharedlibrarybasename := lib$(dirname)
libraryobject         := lib$(dirname).o
sharedlibrary         := lib$(dirname).so
staticlibrary         := lib$(dirname).a


$(depdir): ; @mkdir -p $@

.PHONY: make_dependencies
make_dependencies: $(depdir)
	@for item in $(sources); do $(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $$item -MF .deps/$$item.d; done

$(objects): %.o: %.cpp | $(depdir)
	@echo Compile object and generate dependencies: $@ 
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $*.cpp -MF .deps/$*.d
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c -o $@ 






headerchecks := $(patsubst %.hpp,check-%.hpp,$(headers))

.PHONY: $(headerchecks)
$(headerchecks): check-%.hpp : 
	$(info Check header: $*.hpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.hpp -fsyntax-only

.PHONY: checkheaders
checkheaders: $(headerchecks)




# NOTE: Original recipe for the shared library, now there is just one object
# $(sharedlibrary): $(objects)
# 	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

.all.o: $(sources) .all.cpp | $(depdir)
	@echo Compiling and setting dependencies: $(libraryobject)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM .all.cpp -MF .deps/.all.d
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) .all.cpp -c -o $@ 

$(libraryobject): .all.o
	@echo Create library object: $*
	@cp .all.o $(libraryobject)

$(sharedlibrary): $(libraryobject)
	@echo Shared library: $@
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

$(staticlibrary): $(libraryobject)
	@echo Static library: $*
	@ar rcs $(staticlibrary) $(libraryobject)






-include $(depdir)/.all.d
-include $(dependencies)

*.o .all.o: ./makefile ../makefile ../common.recipe.mk ../common.rules.mk ../common.upkeep.mk


.PHONY: buildobjects buildso builda
buildobjects: $(objects)
buildso:      $(sharedlibrary)
builda:       $(staticlibrary)

.PHONY: build
ifeq ($(OS),Windows_NT)
build: builda
else
build: buildso builda
endif

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

