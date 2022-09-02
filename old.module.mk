SHELL = /bin/bash

# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 



dirname := $(notdir $(CURDIR))
depdir  := .deps

sources      := $(wildcard *.cpp)
headers      := $(wildcard *.hpp)
objects      := $(patsubst %.cpp,%.o,$(sources))
dependencies := $(patsubst %.cpp,$(depdir)/%.d,$(sources))

sharedlibrarybasename := lib$(dirname)
libraryobject         := lib$(dirname).o
sharedlibrary         := lib$(dirname).so
staticlibrary         := lib$(dirname).a



###################################################################################################
# How to compile the .o files and .all.o file 


$(depdir):
	@"mkdir" -p $@
#	@echo $(dirname)

DEPFLAGS = -MT $@ -MMD -MP -MF $(depdir)/$*.d

.PHONY: make_dependencies
make_dependencies: $(depdir)
	@for item in $(sources); do $(CXX) $(CXXFLAGS) $(CPPFLAGS) $$item -MM -MP -MF $(depdir)/$$item.d; done

$(objects): %.o: %.cpp $(depdir)/%.d | $(depdir)
	@echo Compile object and generate dependencies: $@ 
#	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $*.cpp -MF $(depdir)/$*.d #
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS)       $< -c -o $@ $(DEPFLAGS)

.all.o: $(sources) .all.cpp $(depdir)/.all.d | $(depdir)
	@echo Compiling and setting dependencies: $(libraryobject)
#	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM .all.cpp -MF $(depdir)/.all.d
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) .all.cpp -c -o $@  $(DEPFLAGS)

$(depdir)/.all.d:
$(dependencies):

-include $(depdir)/.all.d
-include $(dependencies)


# How to provide the library .o file, and possibly static/shared libraries

$(libraryobject): .all.o
	@echo Library object: $@
	@cp .all.o $(libraryobject)

$(sharedlibrary): $(libraryobject)
	@echo Shared library: $@
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

$(staticlibrary): $(libraryobject)
	@echo Static library: $@
	@ar rcs $(staticlibrary) $(libraryobject)


###################################################################################################
# All object files that are compiled (and their descendants) also depend on the makefiles 

*.o .all.o: ./makefile ../makefile ../common.compile.mk ../common.module.mk ../common.upkeep.mk


###################################################################################################
# Finally, determine the build target depending on whether shared libraries are poossible or not

.PHONY: build buildobjects buildso builda

buildobjects: $(objects)
buildso:      $(sharedlibrary)
builda:       $(staticlibrary)

ifeq ($(OS),Windows_NT)
build: builda
else
build: buildso builda
endif



########################################################################
# List several objects. Read-only

.PHONY:  list_of_objects
.SILENT: list_of_objects
list_of_objects: 
	@echo $(dirname);
	@echo $(staticlibrary);
	@echo $(sharedlibrary);
	@echo $(objects);
	@echo $(dependencies);
	@echo $(build);


########################################################################
# Check whether the header files have correct syntax. Read-only.

headerchecks := $(patsubst %.hpp,check-%.hpp,$(headers))

.PHONY: $(headerchecks) checkheaders

checkheaders: $(headerchecks)

$(headerchecks): check-%.hpp : 
	$(info Check header: $*.hpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.hpp -fsyntax-only









#buildobjects # NOTE: the .o files were required for our .so files originally

# NOTE: Original recipe for the shared library, now there is just one object
# $(sharedlibrary): $(objects)
# 	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

# $(sharedlibrary): #$(objects)
# 	# TODO : create the dependencies of the cpp file 
# 	MYFILE=''; for mycpp in *.cpp; do MYFILE+=$$'#include \"'$$mycpp$$'\"\n'; done; (echo "$$MYFILE") > $(sharedlibrarybasename).code
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -x c++ $(sharedlibrarybasename).code -c -o $(sharedlibrarybasename).o
# 	$(CXX) $(CXXFLAGS) -shared -o $@ $(sharedlibrarybasename).o $(LDLIBS)
