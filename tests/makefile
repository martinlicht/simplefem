SHELL = /bin/sh

.PHONY: default 
default: build

.PHONY: build 
# build is set up below 

projectdir :=.

include ../common.compile.mk 

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
components+=afem
components+=sullivan2D
components+=sullivan3D
components+=whitney2D
components+=whitney3D
components+=solverfem
components+=nullspace


# cleanfiles    :=
# vtkcleanfiles :=
# depcleanfiles :=


include tests.affices.mk

testsdir:=.
projectdir:=../
pathvar:=$(CURDIR)/../

contextdir:=./basic
context:=basic
include tests.rules.mk
#include tests.run.mk
build: basic.tests 

contextdir:=./utility
context:=utility
include tests.rules.mk
#include tests.run.mk
build: utility.tests

contextdir:=./combinatorics
context:=combinatorics
include tests.rules.mk
#include tests.run.mk
build: combinatorics.tests

contextdir:=./operators
context:=operators
include tests.rules.mk
#include tests.run.mk
build: operators.tests

contextdir:=./dense
context:=dense
include tests.rules.mk
#include tests.run.mk
build: dense.tests

contextdir:=./sparse
context:=sparse
include tests.rules.mk
#include tests.run.mk
build: sparse.tests

contextdir:=./solver
context:=solver
include tests.rules.mk
#include tests.run.mk
build: solver.tests

contextdir:=./mesh
context:=mesh
include tests.rules.mk
#include tests.run.mk
build: mesh.tests

contextdir:=./vtk
context:=vtk
include tests.rules.mk
#include tests.run.mk
build: vtk.tests

# contextdir:=./matrixmarket
# context:=matrixmarket
# include tests.rules.mk
# #include tests.run.mk
# build: matrixmarket.tests

contextdir:=./fem
context:=fem
include tests.rules.mk
#include tests.run.mk
build: fem.tests

contextdir:=./solverfem
context:=solverfem
include tests.rules.mk
#include tests.run.mk
build: solverfem.tests

contextdir:=./sullivan2D
context:=sullivan2D
include tests.rules.mk
#include tests.run.mk
build: sullivan2D.tests

contextdir:=./sullivan3D
context:=sullivan3D
include tests.rules.mk
#include tests.run.mk
build: sullivan3D.tests

contextdir:=./whitney2D
context:=whitney2D
include tests.rules.mk
#include tests.run.mk
build: whitney2D.tests

contextdir:=./whitney3D
context:=whitney3D
include tests.rules.mk
#include tests.run.mk
build: whitney3D.tests

contextdir:=./afem
context:=afem
include tests.rules.mk
#include tests.run.mk
build: afem.tests




.PHONY: clean
clean:
	@echo "finished cleaning test directory." 


.PHONY: dependclean
dependclean:
	@echo "finished cleaning test directory dependency information files." 


.PHONY: vtkclean
vtkclean:
	@echo "finished cleaning test directory .vtk files." 






# clean:
# 	echo $(cleanfiles)

# vtkclean:
# 	echo $(vtkcleanfiles)

# depclean:
# 	echo $(depcleanfiles)
	




# clean.components :=$(patsubst %,.clean.%,$(components) )

# .PHONY: $(clean.components)
# $(clean.components): .clean.%: 
# 	@cd ./$* && $(MAKE) --no-print-directory clean

# .PHONY: clean
# clean: $(clean.components)
# 	@$(MAKE) --no-print-directory -f ../common.upkeep.mk clean
# 	@echo "finished cleaning." 


# dependclean.components :=$(patsubst %,.dependclean.%,$(components) )

# .PHONY: $(dependclean.components)
# $(dependclean.components): .dependclean.%: 
# 	@cd ./$* && $(MAKE) --no-print-directory dependclean

# .PHONY: dependclean
# dependclean: $(dependclean.components)
# 	$(MAKE) --no-print-directory -f ../common.upkeep.mk dependclean
# 	@echo "finished cleaning dependency information files." 


# vtkclean.components :=$(patsubst %,.vtkclean.%,$(components) )

# .PHONY: $(vtkclean.components)
# $(vtkclean.components): .vtkclean.%: 
# 	@cd ./$* && $(MAKE) --no-print-directory vtkclean

# .PHONY: vtkclean
# vtkclean: $(vtkclean.components)
# 	$(MAKE) --no-print-directory -f ../common.upkeep.mk vtkclean
# 	@echo "finished cleaning .vtk files." 


# tidy.components :=$(patsubst %,.tidy.%,$(components) )

# .PHONY: $(tidy.components)
# $(tidy.components): .tidy.%: 
# 	@cd ./$* && $(MAKE) --no-print-directory tidy

# .PHONY: tidy
# tidy: $(tidy.components)
# 	$(MAKE) --no-print-directory -f ../common.upkeep.mk tidy


# cppcheck.components :=$(patsubst %,.cppcheck.%,$(components) )

# .PHONY: $(cppcheck.components)
# $(cppcheck.components): .cppcheck.%: 
# 	@cd ./$* && $(MAKE) --no-print-directory cppcheck

# .PHONY: cppcheck
# cppcheck: $(cppcheck.components)
# 	$(MAKE) --no-print-directory -f ../common.upkeep.mk cppcheck


