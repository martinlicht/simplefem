
default: build

include ../../common.recipe.mk 

include ../../common.upkeep.mk

include ../tests.affices.mk

projectdir:=../../
pathvar:=$(shell pwd)/../../

# should be a subdirectory of where the sources are
depdir := .deps

contextdir:=.

context:=$(shell basename $$(pwd))


build: $(context).tests

include ../tests.rules.mk
include ../tests.run.mk

