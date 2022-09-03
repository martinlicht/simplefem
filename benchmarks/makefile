SHELL = /bin/sh

.PHONY:default 
default: build

projectdir :=.

include ../common.compile.mk 

affices:= fem vtk mesh solver dense sparse operators combinatorics utility basic

ifeq ($(OS),Windows_NT)
ending := exe
else
ending := out
endif

sources      := $(wildcard ./*.cpp)
outs         := $(patsubst %.cpp,%.$(ending),$(sources))
runs         := $(patsubst %.cpp,%.run,$(sources))
dependencies := $(patsubst ./%.cpp,./.deps/%.d,$(sources))
depdir       := ./.deps

linkerprefix :=-Wl,
includes := $(patsubst %,-L../%,$(affices))
rpath_t  := $(patsubst %,-rpath=../%,$(affices)) 
rpath    := $(patsubst %,$(linkerprefix)%,$(rpath_t)) 
link     := $(patsubst %,-l%,$(affices))

$(depdir):
	@"mkdir" -p $@

DEPFLAGS = -MT $@ -MF .deps/$*.d -MP -MMD
# We generate dependency files during compilation using the following compiler flags 
# -MT $@ : sets the target of the makefile rule, stripped of any directory components
# -MP    : add dependencies as phony targets
# -MF ...: sets the output file for the rules 
# -MMD   : list headers as a by-product of compiling, excluding system headers

$(outs): ./%.$(ending): ./%.cpp | $(depdir)
	@echo Compiling $@ ...
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes) $(rpath) $(link) -o $@ $(LDLIBS) $(DEPFLAGS)
else
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(includes)          $(link) -o $@ $(LDLIBS) $(DEPFLAGS)
endif
	
-include $(dependencies)

%.$(ending): ./makefile ../makefile ../common.compile.mk

.PHONY: build
build: $(outs)





$(runs): %.run : %.$(ending)
	time ./$<

run: $(runs)




CMD_CLEAN    = rm -f .all.o *.a *.o *.d *.so *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_VTKCLEAN = rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk
CMD_DEPCLEAN = if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean vktclean dependclean

clean: 
	@$(CMD_CLEAN); $(CMD_VTKCLEAN); $(CMD_DEPCLEAN); 
	@echo "Finished cleaning benchmark directory." 

vtkclean: 
	@$(CMD_VTKCLEAN); 
	@echo "Finished cleaning benchmark directory .vtk files." 

dependclean: 
	@$(CMD_DEPCLEAN); 
	@echo "Finished cleaning benchmark directory dependency information files." 


