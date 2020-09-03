SHELL = /bin/sh

default: all

include ../common.recipe.mk 

include ../common.upkeep.mk

include tests.affices.mk


./*/*.out: makefile ../common.recipe.mk


depdir := .deps

projectdir:=../
pathvar:=$(shell pwd)/../



contextdir:=./basic
context:=basic
include testrule.mk

contextdir:=./utility
context:=utility
include testrule.mk

contextdir:=./combinatorics
context:=combinatorics
include testrule.mk

contextdir:=./operators
context:=operators
include testrule.mk

contextdir:=./dense
context:=dense
include testrule.mk

contextdir:=./sparse
context:=sparse
include testrule.mk

contextdir:=./solver
context:=solver
include testrule.mk

contextdir:=./mesh
context:=mesh
include testrule.mk

contextdir:=./vtk
context:=vtk
include testrule.mk

contextdir:=./matrixmarket
context:=matrixmarket
include testrule.mk

contextdir:=./fem
context:=fem
include testrule.mk

contextdir:=./experiments
context:=experiments
include testrule.mk

contextdir:=./afem
context:=afem
include testrule.mk






all: basic.tests 
all: utility.tests
all: combinatorics.tests
all: operators.tests
all: dense.tests
all: sparse.tests
all: solver.tests
all: mesh.tests
all: vtk.tests
all: matrixmarket.tests
all: fem.tests
all: experiments.tests
all: afem.tests



include makefile.runtests


clean: 
	-rm -f *.o *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 
	-rm -f */.deps/*.d
	-rm -f */*.o */*.gch
	-rm -f */*.exe */*.exe.stackdump
	-rm -f */*.out */*.out.stackdump 


