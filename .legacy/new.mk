component must provide:

	projectdir
	componentdir
	include compiler options
	include build recipes
	include clean recipe 
	
project must provide

	projectdir 
	include compiler options 
	
	repeat for every component 
	component dir
	include build recipes 
	include clean recipes 


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
================== Component

default: build

projectdir   :=..
include $(projectdir)/common.compile.mk

componentdir :=.
componentname:=$(shell basename $(CURDIR))

include $(projectdir)/common.recipes.mk



##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
================== Project directory 

projectdir   :=..

include common.compile.mk

componentdir :=./combinatorics
componentname:=combinatorics
include common.recipe.mk

.... repeat ... 







##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################


SHELL = /bin/bash

# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 



$(componentname).depdir  := $(componentdir)/.deps

$(componentname).sources      := $(wildcard $(componentdir)/*.cpp)
$(componentname).headers      := $(wildcard $(componentdir)/*.hpp)
$(componentname).objects      := $(patsubst %.cpp,%.o,$($(componentname).sources))
$(componentname).dependencies := $(patsubst %.cpp,.deps/%.d,$($(componentname).sources))

$(componentname).sharedlibrarybasename := $(componentdir)/lib$(componentname)
$(componentname).libraryobject         := $(componentdir)/lib$(componentname).o
$(componentname).sharedlibrary         := $(componentdir)/lib$(componentname).so
$(componentname).staticlibrary         := $(componentdir)/lib$(componentname).a



###################################################################################################
# How to compile the .o files and .all.o file 


$(componentdir)/$(depdir):
	@"mkdir" -p $@

DEPFLAGS = -MT $@ -MMD -MP -MF .deps/$*.d

.PHONY: make_dependencies $(componentname).make_dependencies
make_dependencies: $(componentname).make_dependencies
make_dependencies: $(componentdir)/$(depdir)
	@for item in $(sources); do $(CXX) $(CXXFLAGS) $(CPPFLAGS) $$item -MM -MP -MF .deps/$$item.d; done

$(componentdir)/$($(componentname).objects): %.o: %.cpp $(componentdir)/$(depdir)/%.d | $(componentdir)/$(depdir)
	@echo Compile object and generate dependencies: $@ 
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS)       $< -c -o $@ $(DEPFLAGS)

$(componentdir)/.all.o: $($(componentname).sources) $(componentdir)/.all.cpp $(componentdir)/$(depdir)/.all.d | $(componentdir)/$(depdir)
	@echo Compiling and setting dependencies: $(libraryobject)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) .all.cpp -c -o $@  $(DEPFLAGS)

$(componentdir)/.deps/.all.d:
$(dependencies):

-include $(depdir)/.all.d
-include $(dependencies)


# How to provide the library .o file, and possibly static/shared libraries

$($(componentname).libraryobject): .all.o
	@echo Library object: $@
	@cp .all.o $(libraryobject)

$($(componentname).sharedlibrary): $($(componentname).libraryobject)
	@echo Shared library: $@
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

$($(componentname).staticlibrary): $($(componentname).libraryobject)
	@echo Static library: $@
	@ar rcs $($(componentname).staticlibrary) $($(componentname).libraryobject)


###################################################################################################
# All object files that are compiled (and their descendants) also depend on the makefiles 

*.o .all.o: ./makefile ../makefile ../common.compile.mk ../common.module.mk ../common.upkeep.mk


###################################################################################################
# Finally, determine the build target depending on whether shared libraries are poossible or not

.PHONY: buildobjects buildso builda
.PHONY: $(componentname).buildobjects $(componentname).buildso $(componentname).builda

$(componentname).buildobjects: $($(componentname).objects)
$(componentname).buildso:      $($(componentname).sharedlibrary)
$(componentname).builda:       $($(componentname).staticlibrary)

buildobjects: $(componentname).buildobjects
buildso:      $(componentname).buildso
builda:       $(componentname).builda


.PHONY: build
.PHONY: $(componentname).build

build: $(componentname).build

ifeq ($(OS),Windows_NT)
$(componentname).build: $(componentname).builda
else
$(componentname).build: $(componentname).builda $(componentname).buildso
endif



########################################################################
# List several objects. Read-only

.PHONY:  list_of_objects $(componentname).list_of_objects
.SILENT: list_of_objects $(componentname).list_of_objects
list_of_objects: $(componentname).list_of_objects
$(componentname).list_of_objects: 
	@echo $(componentdir);
	@echo $($(componentname).staticlibrary);
	@echo $($(componentname).sharedlibrary);
	@echo $($(componentname).objects);
	@echo $($(componentname).dependencies);
	@echo $($(componentname).build);


########################################################################
# Check whether the header files have correct syntax. Read-only.

$(componentname).headerchecks := $(patsubst %.hpp,check-%.hpp,$($(componentname).headers))

.PHONY: checkheaders $($(componentname).headerchecks)

checkheaders: $($(componentname).headerchecks)

$($(componentname).headerchecks): check-%.hpp : 
	$(info Check header: $*.hpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.hpp -fsyntax-only




########################################################################
# Commands for cleaning out numerous files that are not part of the project.
# These remove the following:
# - clean: delete all binary files and temporary output, and the following two.
# - vtkclean: delete all *vtk output files
# - dependclean: delete .deps directories and their content

CMD_CLEAN    := rm -f .all.o *.a *.o *.d *.so *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_VTKCLEAN := rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk
CMD_DEPCLEAN := if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean vktclean dependclean $(component).clean $(component).vktclean $(component).dependclean

clean: $(component).clean
$(component).clean: 
	@-$(CMD_CLEAN); $(CMD_VTKCLEAN); $(CMD_DEPCLEAN); 

vtkclean: $(component).vtkclean
$(component).vtkclean: 
	@-$(CMD_VTKCLEAN); 

dependclean: $(component).dependclean
$(component).dependclean: 
	@-$(CMD_DEPCLEAN); 


########################################################################
# apply clang-tidy to all cpp and hpp files in the directory

.PHONY: tidy $(component).tidy
tidy: $(component).tidy
$(component).tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++2a


########################################################################
# apply cppcheck to all cpp and hpp files in the directory 

.PHONY: cppcheck $(component).cppcheck
cppcheck: $(component).cppcheck
$(component).cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy --enable=warning,style,performance,portability --suppress=duplicateCondition --suppress=assertWithSideEffect --suppress=useStlAlgorithm --std=c++17 -q . ./*pp


########################################################################
# apply cpplint to all cpp and hpp files in the directory 

.PHONY: cpplint $(component).cpplint
cpplint: $(component).cpplint
$(component).cpplint:
	$(error This command is not implemented.)
	( $(projectdir)/Tools/cpplint.py --exclude=$(projectdir)/.private/ --exclude=$(projectdir)/.legacy/ --exclude=$(projectdir)/.playground/ --recursive --filter=-whitespace,-legal,-build/namespace,readability/alt_tokens,readability/todo,readability/inheritance --quiet $(projectdir) ) | sort | uniq -c 2> $(projectdir)/OUTPUT_CPPLINT.txt


########################################################################
# regex several useful things 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues $(component).grepissues
grepissues: $(component).grepissues
$(component).grepissues:
#	@echo Find trailing whitespace...
#	@-grep --line-number --color '\s+$$' -r $(componentdir)/*pp
#	@echo Find non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r $(componentdir)/*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r $(componentdir)/*pp
#	@-grep --line-number --color 'assert(' $(componentdir)/*pp
	@-grep --line-number --color 'cout' $(componentdir)/*pp
#	@-grep --line-number --color -E '\.*[0-9]' $(componentdir)/*pp
	@-grep --line-number --color -E '(0-9)e' $(componentdir)/*pp
	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' $(componentdir)/*pp


