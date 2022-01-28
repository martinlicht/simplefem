SHELL = /bin/sh

default: build

all: build tests benchmarks 
	@echo "Finished all"


# Prints the help message: a list of possible makefile targets

.PHONY: help 
help:
	@echo Use one of the following targets:
	@echo ""
	@echo " help:         Print this help."
	@echo " build:        Build all FEECpp libraries, and compile the tests and benchmarks."
	@echo " tests:        Run all tests for all components."
	@echo " benchmarks:   Perform the benchmarks." 
	@echo " all:          build, test, and benchmark"
	@echo " parameters:   Display the build parameters."
	@echo " checkheaders: Run a static code analysis tool. The particular tool is not specified."
	@echo " check:        Run a static code analysis tool. The particular tool is not specified."
	@echo " clean:        Clean all output from previous builds, tests, and benchmarks,"
	@echo "               including all VTK output"
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
#components+=matrixmarket
components+=fem



# Describe the 'build' target and its subtargets

build.components   :=$(patsubst %,.build.%,   $(components) )
build.components.a :=$(patsubst %,.build.%.a, $(components) )
build.components.so:=$(patsubst %,.build.%.so,$(components) )

.PHONY: build .build.tests .build.benchmarks $(build.components)

build: $(build.components) .build.tests .build.benchmarks

$(build.components): .build.%: 
	@echo Build: $*
	@cd ./$* && $(MAKE) --no-print-directory build

.build.tests:
	@echo Build tests
	@cd ./tests/ && $(MAKE) --no-print-directory build

.build.benchmarks:
	@echo Build Benchmarks
	@cd ./benchmarks/ && $(MAKE) --no-print-directory build

.PHONY: .build.a .build.so $(build.components.a) $(build.components.so) 
.build.a:  $(build.components.a)

.build.so: $(build.components.so)

$(build.components.a): .build.%.a:
	@echo Build A: $*
	@cd ./$* && $(MAKE) --no-print-directory builda

$(build.components.so): .build.%.so: 
	@echo Build SO: $*
	@cd ./$* && $(MAKE) --no-print-directory buildso



# The target 'test' runs all the tests in the test directory 

.PHONY: test
test:
	@cd ./tests && $(MAKE) --no-print-directory run


# The target 'benchmark' runs all the benchmarks in the benchmark directory

.PHONY: benchmarks
benchmarks:
	@cd ./benchmarks && $(MAKE) --no-print-directory 



# Target 'check' is a generic test. Currently, it defaults to 'tidy'

check: tidy


# Check all the headers for stand-alone compilation.
# In particular, check for self-inclusion.

checkheaders.components :=$(patsubst %,.checkheaders.%,$(components) )

.PHONY: $(checkheaders.components)
$(checkheaders.components): .checkheaders.%: 
	@cd ./$* && $(MAKE) --no-print-directory checkheaders

.PHONY: checkheaders
checkheaders: $(checkheaders.components)
	@echo "finished..." 


# Target 'tidy' to call clang-tidy

tidy.components :=$(patsubst %,.tidy.%,$(components) )

.PHONY: $(tidy.components)
$(tidy.components): .tidy.%: 
	@cd ./$* && $(MAKE) --no-print-directory tidy

.PHONY: tidy
tidy: $(tidy.components)
	@$(MAKE) --no-print-directory -f common.upkeep.mk tidy


# Call 'cppcheck'
cppcheck.components :=$(patsubst %,.cppcheck.%,$(components) )

.PHONY: $(cppcheck.components)
$(cppcheck.components): .cppcheck.%: 
	@cd ./$* && $(MAKE) --no-print-directory cppcheck

.PHONY: cppcheck
cppcheck: $(cppcheck.components)
	@$(MAKE) --no-print-directory -f common.upkeep.mk cppcheck

# Call 'cpplint' 

.PHONY: cpplint
cpplint:
	@( ./Tools/cpplint.py --exclude=tests/* --exclude=tests/*/* --exclude=.legacy/* --exclude=.private/* --exclude=.playground/* --recursive --filter=-whitespace,-legal --quiet . ) | sort | uniq -c 2> OUTPUT_CPPLINT.txt

# print the build parameters


.PHONY: parameters 
parameters:
	@make --no-print-directory -f common.recipe.mk parameters
	@true


# Display the value of a variable

print-%:
	$(info [ variable name]: $*)
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info [        origin]: $(origin $*))
	$(info )
	@true

# Display the value of all variables.
.PHONY: printall
printall: $(subst :,\:,$(foreach variable,$(.VARIABLES),print-$(variable)))
#printall: 
#	@$(info $(foreach variable,$(.VARIABLES),print-$(variable)) )
#	@true




# 'Clean' target

clean.components :=$(patsubst %,.clean.%,$(components) )

.PHONY: $(clean.components)
$(clean.components): .clean.%: 
	@cd ./$* && $(MAKE) --no-print-directory clean

.PHONY: clean
clean: $(clean.components)
	@$(MAKE) --no-print-directory -f common.upkeep.mk clean
	@echo "finished cleaning files." 


# Clean all dependency files with this target 

dependclean.components :=$(patsubst %,.dependclean.%,$(components) )

.PHONY: $(dependclean.components)
$(dependclean.components): .dependclean.%: 
	@cd ./$* && $(MAKE) --no-print-directory dependclean

.PHONY: dependclean
dependclean: $(dependclean.components)
	@$(MAKE) --no-print-directory -f common.upkeep.mk dependclean
	@echo "finished cleaning dependency information files." 


# Clean all VTK files

vtkclean.components :=$(patsubst %,.vtkclean.%,$(components) )

.PHONY: $(vtkclean.components)
$(vtkclean.components): .vtkclean.%: 
	@cd ./$* && $(MAKE) --no-print-directory vtkclean

.PHONY: vtkclean
vtkclean: $(vtkclean.components)
	@$(MAKE) --no-print-directory -f common.upkeep.mk vtkclean
	@echo "finished cleaning .vtk files." 




