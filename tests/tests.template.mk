#################################################################################################
# This file is supposed to be included in one of the subdirectories within the 'tests' folder.
# It imports all the necessary targets for building the tests in the respective directory.
#################################################################################################

.PHONY: default
default: build

context    :=$(notdir $(CURDIR))
contextdir :=.
testsdir   :=../
projectdir :=../../
pathvar    :=$(CURDIR)/../../

include ../../common.compile.mk

include ../tests.affices.mk

.PHONY: build
build: $(context).tests

include ../tests.rules.mk

