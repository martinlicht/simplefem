
default: builddefault

include ../../common.recipe.mk 

include ../../common.upkeep.mk

include ../tests.affices.mk

projectdir:=../../
pathvar:=$(shell pwd)/../../

# should be a subdirectory of where the sources are
depdir := .deps

contextdir:=.

context:=$(shell basename $$(pwd))


builddefault: $(context).tests

include ../testrule.mk

