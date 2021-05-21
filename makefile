SHELL = /bin/sh

default: build

all: build tests benchmarks 



# Prints the help message: a list of possible makefile targets

.PHONY: help 
help:
	@echo Use one of the following targets:
	@echo ""
	@echo " help:       Print this help."
	@echo " build:      Build all FEECpp libraries, and compile the tests and benchmarks."
	@echo " tests:      Run all tests for all components."
	@echo " benchmarks: Perform the benchmarks." 
	@echo " all:        build, test, and benchmark"
	@echo " check:      Run a static code analysis tool. The particular tool is not specified."
	@echo " clean:      Clean all output from previous builds, tests, and benchmarks,"
	@echo "             including all VTK output"
	@echo ""
	@echo " The default target is [build]"
	@echo ""
	@echo " build, tests, and benchmark require only a C++17 compiler and GNU Make."







# Describe the different components of the software

components:=
components+=basic
components+=utility
components+=combinatorics
components+=operators
components+=dense
components+=sparse
components+=solver
components+=mesh
components+=vtk
components+=matrixmarket
components+=fem






# Describe the 'build' target and its subtargets

build.components   :=$(patsubst %,.build.%,   $(components) )
build.components.a :=$(patsubst %,.build.%.a, $(components) )
build.components.so:=$(patsubst %,.build.%.so,$(components) )

.PHONY: build .build.tests .build.benchmarks $(build.components)

build: $(build.components) .build.tests .build.benchmarks

$(build.components): .build.%: 
	cd ./$* && $(MAKE) build

.build.tests:
	cd ./tests/ && $(MAKE) build

.build.benchmarks:
	cd ./benchmarks/ && $(MAKE) build

.PHONY: .build.a .build.so $(build.components.a) $(build.components.so) 
.build.a:  $(build.components.a)

.build.so: $(build.components.so)

$(build.components.a): .build.%.a:
	cd ./$* && $(MAKE) builda

$(build.components.so): .build.%.so: 
	cd ./$* && $(MAKE) buildso




# The target 'test' runs all the tests in the test directory 

.PHONY: test
test:
	cd ./tests && $(MAKE) run




# The target 'benchmark' runs all the benchmarks in the benchmark directory

.PHONY: benchmarks
benchmarks:
# 	cd ./benchmarks && $(MAKE) 





# Target 'check' is a generic test. Currently, it defaults to 'tidy'

check: tidy





# 'Clean' target

clean.components :=$(patsubst %,.clean.%,$(components) )

.PHONY: $(clean.components)
$(clean.components): .clean.%: 
	cd ./$* && $(MAKE) clean

.PHONY: clean
clean: $(clean.components)
	$(MAKE) -f common.upkeep.mk clean
	@echo "finished cleaning .vtk files." 



# Clean all dependency files with this target 

dependclean.components :=$(patsubst %,.dependclean.%,$(components) )

.PHONY: $(dependclean.components)
$(dependclean.components): .dependclean.%: 
	cd ./$* && $(MAKE) dependclean

.PHONY: dependclean
dependclean: $(dependclean.components)
	$(MAKE) -f common.upkeep.mk dependclean
	@echo "finished cleaning dependency information files." 







# Clean all VTK files

vtkclean.components :=$(patsubst %,.vtkclean.%,$(components) )

.PHONY: $(vtkclean.components)
$(vtkclean.components): .vtkclean.%: 
	cd ./$* && $(MAKE) vtkclean

.PHONY: vtkclean
vtkclean: $(vtkclean.components)
	$(MAKE) -f common.upkeep.mk vtkclean
	@echo "finished cleaning." 




# Target 'tidy' to call clang-tidy

tidy.components :=$(patsubst %,.tidy.%,$(components) )

.PHONY: $(tidy.components)
$(tidy.components): .tidy.%: 
	cd ./$* && $(MAKE) tidy

.PHONY: tidy
tidy: $(tidy.components)
	$(MAKE) -f common.upkeep.mk tidy



# Call 'cppcheck'
cppcheck.components :=$(patsubst %,.cppcheck.%,$(components) )

.PHONY: $(cppcheck.components)
$(cppcheck.components): .cppcheck.%: 
	cd ./$* && $(MAKE) cppcheck

.PHONY: cppcheck
cppcheck: $(cppcheck.components)
	$(MAKE) -f common.upkeep.mk cppcheck




# Call 'cpplint' 

.PHONY: cpplint
cpplint:
	( ./Tools/cpplint.py --exclude=tests/* --exclude=tests/*/* --exclude=.legacy/* --exclude=.private/* --exclude=.playground/* --recursive --filter=-whitespace,-legal --quiet . ) | sort | uniq -c > OUTPUT_CPPLINT.txt






# .PHONY: tidy
# tidy: 
# 	cd ./utility && $(MAKE) tidy
# 	cd ./combinatorics && $(MAKE) tidy
# 	cd ./operators && $(MAKE) tidy
# 	cd ./dense && $(MAKE) tidy
# 	cd ./sparse && $(MAKE) tidy
# 	cd ./solver && $(MAKE) tidy
# 	cd ./mesh && $(MAKE) tidy
# 	cd ./vtk && $(MAKE) tidy
# 	cd ./matrixmarket && $(MAKE) tidy
# 	cd ./fem && $(MAKE) tidy
# 	cd ./tests && $(MAKE) tidy
# 	$(MAKE) -f makefile.tidy tidy
# 	@echo "finished tidying." 

# 
# .PHONY: check
# check: 
# 	cd ./utility && $(MAKE) check
# 	cd ./combinatorics && $(MAKE) check
# 	cd ./operators && $(MAKE) check
# 	cd ./dense && $(MAKE) check
# 	cd ./sparse && $(MAKE) check
# 	cd ./solver && $(MAKE) check
# 	cd ./mesh && $(MAKE) check
# 	cd ./vtk && $(MAKE) check
# 	cd ./matrixmarket && $(MAKE) check
# 	cd ./fem && $(MAKE) check
# # 	cd ./tests && $(MAKE) check
# 	$(MAKE) -f makefile.check check
# 	@echo "finished checking." 









# .PHONY: check
# check:
# 	cppcheck $(CHECK_OPTION) . -ilegacy/ -iplayground/
# 	cppcheck $(CHECK_OPTION) basic.hpp
# 	cppcheck $(CHECK_OPTION) basic/*.?pp
# 	cppcheck $(CHECK_OPTION) combinatorics/*.?pp
# 	cppcheck $(CHECK_OPTION) operators/*.?pp
# 	cppcheck $(CHECK_OPTION) dense/*.?pp
# 	cppcheck $(CHECK_OPTION) sparse/*.?pp
# 	cppcheck $(CHECK_OPTION) solver/*.?pp
# 	cppcheck $(CHECK_OPTION) mesh/*.?pp
# 	cppcheck $(CHECK_OPTION) */*.?pp



# 	cd ./utility && $(MAKE) cpplint
# 	cd ./combinatorics && $(MAKE) cpplint
# 	cd ./operators && $(MAKE) cpplint
# 	cd ./dense && $(MAKE) cpplint
# 	cd ./sparse && $(MAKE) cpplint
# 	cd ./solver && $(MAKE) cpplint
# 	cd ./mesh && $(MAKE) cpplint
# 	cd ./vtk && $(MAKE) cpplint
# 	cd ./matrixmarket && $(MAKE) cpplint
# 	cd ./fem && $(MAKE) cpplint
# # 	cd ./tests && $(MAKE) cpplint
# 	$(MAKE) -f makefile.cpplint cpplint
# 	@echo "finished cpplinting." 






# .PHONY: foobuild 
# foobuild: 
# 	for i in $(components); do if [ -d $$i ]; then cd $$i; $(MAKE) build; cd ./..; fi; done

# .PHONY: fooclean 
# fooclean: 
# 	for i in $(components); do if [ -d $$i ]; then cd $$i; $(MAKE) clean; cd ./..; fi; done


# 
# buildobjects:
# 	cd ./utility && $(MAKE) buildobjects
# 	cd ./combinatorics && $(MAKE) buildobjects 
# 	cd ./operators && $(MAKE) buildobjects
# 	cd ./dense && $(MAKE) buildobjects
# 	cd ./sparse && $(MAKE) buildobjects
# 	cd ./solver && $(MAKE) buildobjects
# 	cd ./mesh && $(MAKE) buildobjects
# 	cd ./vtk && $(MAKE) buildobjects
# 	cd ./matrixmarket && $(MAKE) buildobjects
# 	# cd ./fem && $(MAKE)
# 	echo "finished building objects" 
# 
# buildso:
# 	cd ./utility && $(MAKE) buildso
# 	cd ./combinatorics && $(MAKE) buildso
# 	cd ./operators && $(MAKE) buildso
# 	cd ./dense && $(MAKE) buildso
# 	cd ./sparse && $(MAKE) buildso
# 	cd ./solver && $(MAKE) buildso
# 	cd ./mesh && $(MAKE) buildso
# 	cd ./vtk && $(MAKE) buildso
# 	cd ./matrixmarket && $(MAKE) buildso
# 	# cd ./fem && $(MAKE)
# 	echo "finished building shared libraries." 

# 
# lib:
# 	cd ./utility && $(MAKE) 
# 	cd ./combinatorics && $(MAKE) 
# 	cd ./operators && $(MAKE) 
# 	cd ./dense && $(MAKE) 
# 	cd ./sparse && $(MAKE) 
# 	cd ./solver && $(MAKE) 
# 	cd ./mesh && $(MAKE) 
# 	cd ./vtk && $(MAKE)
# 	cd ./matrixmarket && $(MAKE)
# 	# cd ./fem && $(MAKE)
# 	echo "finished building objects and shared libraries." 

# 
# .PHONY: test 
# test:
# 	cd ./tests && $(MAKE)
# 	echo "finished building tests."

# 
# .PHONY: clean
# clean: 
# 	cd ./utility && $(MAKE) clean
# 	cd ./combinatorics && $(MAKE) clean
# 	cd ./operators && $(MAKE) clean
# 	cd ./dense && $(MAKE) clean
# 	cd ./sparse && $(MAKE) clean
# 	cd ./solver && $(MAKE) clean
# 	cd ./mesh && $(MAKE) clean
# 	cd ./vtk && $(MAKE) clean
# 	cd ./matrixmarket && $(MAKE) clean
# 	cd ./fem && $(MAKE) clean
# 	cd ./tests && $(MAKE) clean
# 	$(MAKE) -f common.upkeep.mk clean
# 	@echo "finished cleaning." 
# 
# 
# .PHONY: dependclean
# dependclean: 
# 	cd ./utility && $(MAKE) dependclean
# 	cd ./combinatorics && $(MAKE) dependclean
# 	cd ./operators && $(MAKE) dependclean
# 	cd ./dense && $(MAKE) dependclean
# 	cd ./sparse && $(MAKE) dependclean
# 	cd ./solver && $(MAKE) dependclean
# 	cd ./mesh && $(MAKE) dependclean
# 	cd ./vtk && $(MAKE) dependclean
# 	cd ./matrixmarket && $(MAKE) dependclean
# 	cd ./fem && $(MAKE) dependclean
# 	cd ./tests && $(MAKE) dependclean
# 	$(MAKE) -f common.upkeep.mk dependclean
# 	@echo "finished cleaning dependencies." 
