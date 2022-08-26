##########################################################################
##########################################################################
##########################################################################
################## Project directory makefile 
SHELL = /bin/sh


# The default target builds all files  

.PHONY: default all build tests benchmarks 

default: build

all: build tests benchmarks 
	@echo "Finished all"



# Prints the help message: a list of possible makefile targets

.PHONY: help 
help:
	@echo "Use one of the following targets:"
	@echo ""
	@echo " help:         Print this help."
	@echo " build:        Build all FEECpp libraries, and compile the tests and benchmarks."
	@echo " tests:        Run all tests for all modules."
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
	@echo " build, tests, and benchmark require only a C++14 compiler and GNU Make."










# Describe the 'build' target and its subtargets
# Recap all this .... 

build.modules   :=$(patsubst %,.build.%,   $(modules) )
build.modules.a :=$(patsubst %,.build.%.a, $(modules) )
build.modules.so:=$(patsubst %,.build.%.so,$(modules) )

build: .build.modules .build.tests .build.benchmarks

.build.modules: $(build.modules)

.build.tests:
	@echo Build tests
	@cd ./tests/ && $(MAKE) --no-print-directory build

.build.benchmarks:
	@echo Build Benchmarks
	@cd ./benchmarks/ && $(MAKE) --no-print-directory build

.build.a:  $(build.modules.a)

.build.so: $(build.modules.so)

$(build.modules): .build.%: %
	@echo Build: $*
	@cd ./$* && $(MAKE) --no-print-directory build

$(build.modules.a): .build.%.a:
	@echo Build A: $*
	@cd ./$* && $(MAKE) --no-print-directory builda

$(build.modules.so): .build.%.so: 
	@echo Build SO: $*
	@cd ./$* && $(MAKE) --no-print-directory buildso

.PHONY: build .build.tests .build.benchmarks $(build.modules)
.PHONY: .build.a .build.so $(build.modules.a) $(build.modules.so) 



# # The target 'test' runs all the tests in the test directory 

# test:
# 	@cd ./tests && $(MAKE) --no-print-directory run

# .PHONY: test


# # The target 'benchmark' runs all the benchmarks in the benchmark directory

# benchmarks:
# 	@cd ./benchmarks && $(MAKE) --no-print-directory 

# .PHONY: benchmarks









# Describe the different modules of the software
# Provide all targets for the different modules in building the modules 


modules:=
modules+=basic
modules+=utility
modules+=combinatorics
modules+=operators
modules+=dense
modules+=sparse
modules+=solver
modules+=mesh
modules+=vtk
#modules+=matrixmarket
modules+=fem


projectdir   :=.

include common.compile.mk

moddir :=./basic
module:=basic
include common.module.mk

moddir :=./utility
module:=utility
include common.module.mk

moddir :=./combinatorics
module:=combinatorics
include common.module.mk

moddir :=./operators
module:=operators
include common.module.mk

moddir :=./dense
module:=dense
include common.module.mk

moddir :=./sparse
module:=sparse
include common.module.mk

moddir :=./solver
module:=solver
include common.module.mk

moddir :=./mesh
module:=mesh
include common.module.mk

moddir :=./vtk
module:=vtk
include common.module.mk

moddir :=./fem
module:=fem
include common.module.mk






########################################################################
# Apply cpplint to all cpp and hpp files in the directory. Read-only. 

.PHONY: cpplint $(module).cpplint
cpplint: $(module).cpplint
$(module).cpplint:
	$(error This command is not implemented.)
	( $(projectdir)/Tools/cpplint.py --exclude=$(projectdir)/.private/ --exclude=$(projectdir)/.legacy/ --exclude=$(projectdir)/.playground/ --recursive --filter=-whitespace,-legal,-build/namespace,readability/alt_tokens,readability/todo,readability/inheritance --quiet $(projectdir) ) | sort | uniq -c 2> $(projectdir)/OUTPUT_CPPLINT.txt



# print the build parameters

.PHONY: parameters 
parameters:
	@$(MAKE) --no-print-directory -f common.compile.mk parameters
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



