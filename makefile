##########################################################################
##########################################################################
##########################################################################
################## Project directory makefile ############################
##########################################################################
##########################################################################
##########################################################################

SHELL = /bin/sh

# TODO: all phony targets should be declared as such in all files 

# The main targets are 
#  - all 
#  - build 
#  - tests
# The default target builds all files  

default: build

all: build tests
	@echo "Finished all"

.PHONY: default all build tests




# Prints the help message: a list of possible makefile targets

.PHONY: help 
help:
	@echo "Use one of the following targets:"
	@echo ""
	@echo " help:         Print this help."
	@echo " build:        Build all FEECpp libraries, and compile the tests."
	@echo " tests:        Run all tests for all modules."
	@echo " all:          build and tests"
	@echo " parameters:   Display the build parameters."
	@echo " checkheaders: Run a static code analysis tool. The particular tool is not specified."
	@echo " checksources: Run a static code analysis tool. The particular tool is not specified."
	@echo " check:        Run a static code analysis tool. The particular tool is not specified."
	@echo " clean:        Clean all output from previous builds, tests,"
	@echo "               including all VTK output"
	@echo ""
	@echo " The default target is [build]"
	@echo ""
	@echo " build and tests require a C++14 compiler and GNU Make."











################################################################## 
################################################################## 
################################################################## 
################################################################## 

# Describe the different modules of the software
# Define all targets for the different modules in building the modules 

modules:=
modules+=external
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

.PHONY: .buildmodules
.buildmodules: $(patsubst %,%.build,$(modules))

projectdir :=.

include common.compile.mk

moddir :=./external
module:=external
include common.module.mk

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



################################################################## 
################################################################## 
################################################################## 
################################################################## 

# The 'build' target depends on the builds of modules, tests, and benchmarks


.PHONY: build .buildmodules .buildtests # .buildbenchmarks

build: .buildmodules .buildtests # .buildbenchmarks

.buildmodules:
	@echo Built modules 
	
.buildtests: .buildmodules
	@cd ./tests/ && $(MAKE) --no-print-directory build
	@echo Built tests





# The target 'test' runs all the tests in the test directory 

test:
	@cd ./tests && $(MAKE) --no-print-directory run
	@echo Run all tests

# benchmarks:
# 	@cd ./benchmarks && $(MAKE) --no-print-directory 

.PHONY: test # benchmarks






# Upkeep targets that remove clutter and temporary files 
# Also commands for source code checking 
# More targets are defined in each module 

.PHONY: clean
clean: 
	@cd ./tests && $(MAKE) --no-print-directory clean
	@rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk ./*.svg ./*/*.svg ./*/*/*.svg ./*.tex ./*/*.tex ./*/*/*.tex
	@echo "Finished cleaning."

.PHONY: outputclean
outputclean:
	@cd ./tests && $(MAKE) --no-print-directory outputclean
	@rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk ./*.svg ./*/*.svg ./*/*/*.svg ./*.tex ./*/*.tex ./*/*/*.tex
	@echo "Finished cleaning output files."

.PHONY: dependclean
dependclean:
	@cd ./tests && $(MAKE) --no-print-directory dependclean
	@echo "Finished cleaning dependency information files."

.PHONY: tidy
tidy:
	@cd ./tests && $(MAKE) --no-print-directory tidy

.PHONY: cppcheck
cppcheck:
	@cd ./tests && $(MAKE) --no-print-directory cppcheck

.PHONY: check
check:
	@cd ./tests && $(MAKE) --no-print-directory check




###############################################################################################
####   Apply cpplint to all cpp and hpp files in the entire source directory. Read-only.   ####
###############################################################################################

.PHONY: cpplint
cpplint:
	( $(projectdir)/.Tools/cpplint.py \
	--exclude=$(projectdir)/PACKAGE/* --exclude=$(projectdir)/.private/* --exclude=$(projectdir)/.legacy/* --exclude=$(projectdir)/.stuff/* --exclude=$(projectdir)/.CUDA/* --exclude=$(projectdir)/.playground/* --exclude=$(projectdir)/external/* \
	--filter=-runtime/int,-runtime/references,-whitespace,-legal,-build/namespace,-readability/alt_tokens,-readability/fn_size,-readability/todo,-readability/inheritance,-readability/braces,-runtime/arrays,-build/header_guard,-build/include,-build/c++11 \
	--recursive --quiet $(projectdir) 2>&1 ) | \
	sort | uniq | \
	cat > $(projectdir)/OUTPUT_CPPLINT.txt; \
	cat $(projectdir)/OUTPUT_CPPLINT.txt




###############################################################################################
#######    Display the value of all variables defined after parsing this makefile     #########
###############################################################################################

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
