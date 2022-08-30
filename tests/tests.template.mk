# 
# This file is supposed to be included in one of the subdirectories 
# within the 'tests' folder. It imports all the stuff it needs.
# 

PHONY: default 
default: build

include ../../common.compile.mk 

include ../tests.affices.mk

projectdir:=../../
pathvar:=$(CURDIR)/../../

# should be a subdirectory of where the sources are
# TODO: remove depdir from the discussion 
depdir := .deps


context    :=$(notdir $(CURDIR))
contextdir :=.

PHONY: build 
build: $(context).tests

include ../tests.rules.mk

# clean:
# 	echo $(cleanfiles)

# vtkclean:
# 	echo $(vtkcleanfiles)

# depclean:
# 	echo $(depcleanfiles)