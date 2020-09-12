SHELL = /bin/sh

default: all

include ../common.recipe.mk 

include ../common.upkeep.mk

include tests.affices.mk


%.out: makefile ../common.recipe.mk


depdir := .deps

projectdir:=../
pathvar:=$(shell pwd)/../



contextdir:=./basic
context:=basic
include tests.rules.mk
include tests.run.mk
all: basic.tests 

contextdir:=./utility
context:=utility
include tests.rules.mk
include tests.run.mk
all: utility.tests

contextdir:=./combinatorics
context:=combinatorics
include tests.rules.mk
include tests.run.mk
all: combinatorics.tests

contextdir:=./operators
context:=operators
include tests.rules.mk
include tests.run.mk
all: operators.tests

contextdir:=./dense
context:=dense
include tests.rules.mk
include tests.run.mk
all: dense.tests

contextdir:=./sparse
context:=sparse
include tests.rules.mk
include tests.run.mk
all: sparse.tests

contextdir:=./solver
context:=solver
include tests.rules.mk
include tests.run.mk
all: solver.tests

contextdir:=./mesh
context:=mesh
include tests.rules.mk
include tests.run.mk
all: mesh.tests

contextdir:=./vtk
context:=vtk
include tests.rules.mk
include tests.run.mk
all: vtk.tests

contextdir:=./matrixmarket
context:=matrixmarket
include tests.rules.mk
include tests.run.mk
all: matrixmarket.tests

contextdir:=./fem
context:=fem
include tests.rules.mk
include tests.run.mk
all: fem.tests

contextdir:=./experiments
context:=experiments
include tests.rules.mk
include tests.run.mk
all: experiments.tests

contextdir:=./afem
context:=afem
include tests.rules.mk
include tests.run.mk
all: afem.tests




# clean: 
# 	-rm -f *.o *.gch
# 	-rm -f *.exe *.exe.stackdump
# 	-rm -f *.out *.out.stackdump 
# 	-rm -f */.deps/*.d
# 	-rm -f */*.o */*.gch
# 	-rm -f */*.exe */*.exe.stackdump
# 	-rm -f */*.out */*.out.stackdump 


