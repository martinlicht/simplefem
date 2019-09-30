
SHELL = /bin/sh

default: all

dirname := $(notdir $(shell pwd))
depdir := .deps

$(depdir): ; @mkdir -p $@

sources := $(wildcard *.cpp)
objects := $(patsubst %.cpp,%.o,$(sources))
dependencies := $(patsubst %.cpp,.deps/%.d,$(sources))

library := lib$(dirname).so

all: $(objects) $(library)

$(library): $(objects)
	$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)


# $(dependencies): .deps/%.d: %.cpp | $(depdir)
# g++ -MM $*.cpp -MF .deps/$*.d

$(objects): %.o: %.cpp | $(depdir)
	@g++ -MM $*.cpp -MF .deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c -o $@ 

.PHONY:  list_of_objects
.SILENT: list_of_objects
list_of_objects: 
		@echo $(dirname);
		@echo $(library);
		@echo $(objects);
		@echo $(dependencies);


-include $(dependencies)
# -include $(sources:.cpp=.d) .deps/

*.o: ../common.recipe.mk ../common.rules.mk ../common.upkeep.mk ./makefile




