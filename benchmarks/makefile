SHELL = /bin/sh

default: build

include ../common.compile.mk 
include ../common.upkeep.mk

affices:= fem vtk mesh solver dense sparse operators combinatorics basic

depdir := .deps


sources      := $(wildcard ./*.cpp)
outs         := $(patsubst %.cpp,%.out,$(sources))
runs         := $(patsubst %.cpp,%.run,$(sources))
dependencies := $(patsubst ./%.cpp,./.deps/%.d,$(sources))
depdir       := ./$(depdir)

linkerprefix :=-Wl,
includes := $(patsubst %,-L../%,$(affices))
rpath_t  := $(patsubst %,-rpath=../%,$(affices)) 
rpath    := $(patsubst %,$(linkerprefix)%,$(rpath_t)) 
link     := $(patsubst %,-l%,$(affices))

$(depdir): ; @mkdir -p $@

$(outs): ./%.out: ./%.cpp | $(depdir)
	@echo Compiling $@ ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM ./$*.cpp -MT $@ -MF $(depdir)/$*.d
ifeq ($(OS),Windows_NT)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes) $(rpath) $(link) -o $@ $(LDLIBS)
else
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes) $(rpath) $(link) -o $@ $(LDLIBS)
endif
	
-include $(dependencies)

%.out: ./makefile ../makefile ../common.compile.mk ../common.upkeep.mk

build: $(outs)





$(runs): %.run : %.out
	time ./$<

run: $(runs)




