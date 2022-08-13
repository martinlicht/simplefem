SHELL = /bin/sh

default: build

include ../common.compile.mk 
include ../common.upkeep.mk

affices:= fem vtk mesh solver dense sparse operators combinatorics basic

depdir := .deps

ifeq ($(OS),Windows_NT)
ending := exe
else
ending := out
endif

sources      := $(wildcard ./*.cpp)
outs         := $(patsubst %.cpp,%.$(ending),$(sources))
runs         := $(patsubst %.cpp,%.run,$(sources))
dependencies := $(patsubst ./%.cpp,./.deps/%.d,$(sources))
depdir       := ./$(depdir)

linkerprefix :=-Wl,
includes := $(patsubst %,-L../%,$(affices))
rpath_t  := $(patsubst %,-rpath=../%,$(affices)) 
rpath    := $(patsubst %,$(linkerprefix)%,$(rpath_t)) 
link     := $(patsubst %,-l%,$(affices))

$(depdir):
	@mkdir -p $@

$(outs): ./%.$(ending): ./%.cpp | $(depdir)
	@echo Compiling $@ ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM ./$*.cpp -MT $@ -MF $(depdir)/$*.d
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes) $(rpath) $(link) -o $@ $(LDLIBS)
else
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes)          $(link) -o $@ $(LDLIBS)
endif
	
-include $(dependencies)

%.$(ending): ./makefile ../makefile ../common.compile.mk ../common.upkeep.mk

build: $(outs)





$(runs): %.run : %.$(ending)
	time ./$<

run: $(runs)




