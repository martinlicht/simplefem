
SHELL = /bin/sh


# The default target builds all files  

default: build

all: build tests benchmarks 
	@echo "Finished all"

.PHONY: default all 


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
	@echo " build, tests, and benchmark require only a C++20 compiler and GNU Make."


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

build: .build.modules .build.tests .build.benchmarks

.build.modules: $(build.components)

.build.tests:
	@echo Build tests
	@cd ./tests/ && $(MAKE) --no-print-directory build

.build.benchmarks:
	@echo Build Benchmarks
	@cd ./benchmarks/ && $(MAKE) --no-print-directory build

.build.a:  $(build.components.a)

.build.so: $(build.components.so)

$(build.components): .build.%: %
	@echo Build: $*
	@cd ./$* && $(MAKE) --no-print-directory build

$(build.components.a): .build.%.a:
	@echo Build A: $*
	@cd ./$* && $(MAKE) --no-print-directory builda

$(build.components.so): .build.%.so: 
	@echo Build SO: $*
	@cd ./$* && $(MAKE) --no-print-directory buildso

.PHONY: build .build.tests .build.benchmarks $(build.components)
.PHONY: .build.a .build.so $(build.components.a) $(build.components.so) 



# The target 'test' runs all the tests in the test directory 

test:
	@cd ./tests && $(MAKE) --no-print-directory run

.PHONY: test


# The target 'benchmark' runs all the benchmarks in the benchmark directory

benchmarks:
	@cd ./benchmarks && $(MAKE) --no-print-directory 

.PHONY: benchmarks




# Call 'cpplint' on the entire source stack. Read-only 

cpplint:
#	@ cpplint --exclude=tests/* --exclude=tests/*/* --exclude=.legacy/* --exclude=.private/* --exclude=.playground/* --exclude=tests/.legacy/ --recursive --filter=-readability/todo,-build/header_guard,-build/include,-readability/alt_tokens,-whitespace,-legal --quiet . 2>&1 | sort | uniq -c > OUTPUT_CPPLINT.txt ; cat OUTPUT_CPPLINT.txt
	@ ./Tools/cpplint.py --exclude=tests/* --exclude=tests/*/* --exclude=.legacy/* --exclude=.private/* --exclude=.playground/* --exclude=tests/.legacy/ --recursive --filter=-readability/todo,-build/header_guard,-build/include,-readability/alt_tokens,-whitespace,-legal --quiet . 2>&1 | sort | uniq -c > OUTPUT_CPPLINT.txt ; cat OUTPUT_CPPLINT.txt

.PHONY: cpplint


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



