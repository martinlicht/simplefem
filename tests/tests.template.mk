# 
# This file is supposed to be included in one of the subdirectories 
# within the 'tests' folder. It imports all the stuff it needs.
# 


default: build

include ../../common.compile.mk 

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
# include ../tests.run.mk

