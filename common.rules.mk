
SHELL = /bin/sh

default: all

dirname := $(notdir $(shell pwd))
depdir := .deps

$(depdir): ; @mkdir -p $@

sources := $(wildcard *.cpp)
objects := $(patsubst %.cpp,%.o,$(sources))
dependencies := $(patsubst %.cpp,.deps/%.d,$(sources))

sharedlibrary := lib$(dirname).so

buildobjects: $(objects)
buildso: $(sharedlibrary)

# all: $(objects) $(sharedlibrary)
all: buildobjects buildso

$(sharedlibrary): $(objects)
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)


.PHONY: make_dependencies
make_dependencies:
	for item in $(sources); do g++ -MM $$item -MF .deps/$*.d; done

$(objects): %.o: %.cpp | $(depdir)
	@g++ -MM $*.cpp -MF .deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c -o $@ 

.PHONY:  list_of_objects
.SILENT: list_of_objects
list_of_objects: 
		@echo $(dirname);
		@echo $(sharedlibrary);
		@echo $(objects);
		@echo $(dependencies);


-include $(dependencies)
# -include $(sources:.cpp=.d) .deps/

*.o: ../common.recipe.mk ../common.rules.mk ../common.upkeep.mk ./makefile




