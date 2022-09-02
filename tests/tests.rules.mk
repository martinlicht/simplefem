# # Example usage:
# # 
# # 
# # include ../../common.compile.mk 
# # 
# # include ../../common.upkeep.mk
# # 
# # depdir := .deps
# # 
# # context=dense
# # contextdir=.
# # affix.$(context)=dense sparse operators combinatorics
# # projectdir=../../
# # pathvar=$(shell pwd)/../../


################################################################################
# EXPECTED VARIABLES:
# - projectdir : path/to/project/directory
# - context    : path/to/test/component/directory
# - contextdir : name of the test component
ifndef projectdir
$(error Expect 'projectdir')
endif
ifndef context
$(error Expect 'context')
endif
ifndef contextdir
$(error Expect 'contextdir')
endif

$(context).%: mycontext    := $(context)
$(context).%: mycontextdir := $(contextdir)

######################################################################################
# Decide whether to use .exe or .out as executable file ending 

ifeq ($(OS),Windows_NT)
ending := exe
else
ending := out
endif



######################################################################################
# Set the variables for this file 
# determine whether to use static or dynamic linking 

$(context).sources := $(sort $(wildcard $(contextdir)/*.cpp))

$(context).outs    := $(patsubst %.cpp,%.$(ending),$($(context).sources))

$(context).depdir  := $(contextdir)/.deps

$(context).dependencies := $(patsubst $(contextdir)/%.cpp,$(contextdir)/.deps/%.d,$($(context).sources))

$(context).include := $(patsubst %,-L$(projectdir)/%,$(affix.$(context)))
linkerprefix       :=-Wl,
$(context).rpath_t := $(patsubst %,-rpath=$(pathvar)/%,$(affix.$(context))) 
$(context).rpath   := $(patsubst %,$(linkerprefix)%,$($(context).rpath_t)) 
$(context).lib     := $(patsubst %,-l%,$(affix.$(context)))
$(context).alib    := $(patsubst %,-l:lib%.a,$(affix.$(context)))
$(context).solib   := $(patsubst %,-l:lib%.so,$(affix.$(context)))

ifeq ($(LINKINGTYPE),static)
$(context).mylib := $($(context).alib)
else ifeq ($(LINKINGTYPE),dynamic)
$(context).mylib := $($(context).solib)
else ifeq ($(LINKINGTYPE),unspecified)
$(context).mylib := $($(context).lib)
else
$(error No linking mode recognized: $(LINKINGTYPE)) 
endif


###################################################################################################
# All executables compiled depend on the makefiles 

%.$(ending): $(testsdir)/makefile $(testsdir)/tests.rules.mk $(testsdir)/tests.affices.mk $(projectsdir)/common.compile.mk

##########################################################################
# Specific definitions and recipes on how to compile the executables 

.PHONY: $($(context).depdir)
$($(context).depdir):
	@"mkdir" -p $@

DEPFLAGS = -MT $@ -MF $($(mycontext).depdir)/$*.d -MP -MMD
# We generate dependency files during compilation using the following compiler flags 
# -MT $@ : sets the target of the makefile rule, stripped of any directory components
# -MP    : add dependencies as phony targets
# -MF ...: sets the output file for the rules 
# -MMD   : list headers as a by-product of compiling, excluding system headers


$($(context).outs): mycontext    := $(context)
$($(context).outs): mycontextdir := $(contextdir)
$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/%.cpp | $($(context).depdir)
#	@ echo link type:   $(LINKINGTYPE)
#	@ echo context:     $(mycontext)
#	@ echo context dir: $(mycontextdir)
#	@ echo target:      $@
#	@ echo target:      $<
#	@ echo target:      $^
#	@ echo contextdir:  $(mycontextdir)
#	@ echo depdir:      $($(mycontext).depdir)
#	@ echo include:     $($(mycontext).include)
#	@ echo rpath:       $($(mycontext).rpath)
#	@ echo lib:         $($(mycontext).lib)
	@echo Compiling $@ ...
#	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -std=c++2a -MT $@ -MF $($(mycontext).depdir)/$*.d -MM $(mycontextdir)/$*.cpp
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include) $($(mycontext).rpath) $($(mycontext).mylib) -o $@ $(LDLIBS) $(DEPFLAGS)
else
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include)                       $($(mycontext).mylib) -o $@ $(LDLIBS) $(DEPFLAGS)
endif
	@touch $@
	
$($(context).dependencies):

-include $($(context).dependencies)

$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/*.mk $(projectdir)/tests/*.mk

.PHONY: $(context).tests
$(context).tests: $($(context).outs)





###################################################
# How to automatically run the executables

$(context).runs        := $(patsubst %.cpp,%.run,$($(context).sources))

run: $(context).run
$(context).run: $($(context).runs)
$($(context).runs): %.run : %.$(ending)
	./$< 

$(context).silent_runs := $(patsubst %.cpp,%.silent_run,$($(context).sources))

silent_run: $(context).silent_run
$(context).silent_run: $($(context).silent_runs)
$($(context).silent_runs): %.silent_run : %.$(ending)
	./$< > /dev/null 

.PHONY: run silent_run $(context).run $(context).silent_run $($(context).runs) $($(context).silent_runs)

# # 2> /dev/null






# TODO: clean out the stuff below and adapt to the test directory


########################################################################
# Apply clang-tidy to all cpp and hpp files in the directory. Read-only.

.PHONY: tidy $(context).tidy
tidy: $(context).tidy
$(context).tidy:
	clang-tidy $(mycontextdir)/*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++2a


########################################################################
# Apply cppcheck to all cpp and hpp files in the directory. Read-only. 

.PHONY: cppcheck $(context).cppcheck
cppcheck: $(context).cppcheck
$(context).cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy \
	--enable=warning,style,performance,portability --suppress=duplicateCondition \
	--suppress=assertWithSideEffect --suppress=useStlAlgorithm \
	--std=c++17 -q $(mycontextdir)/*pp


########################################################################
# Regex several useful things. Read-only. 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues $(context).grepissues
grepissues: $(context).grepissues
# $(context).grepissues: mycontext := $(context)
# $(context).grepissues: mycontextdir := $(contextdir)
$(context).grepissues:
#	@echo Search trailing whitespace...
#	@-grep --line-number --color '\s+$$' -r $(mycontextdir)/*pp
#	@echo Search non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r $(mycontextdir)/*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r $(mycontextdir)/*pp
#	@echo Find standard asserts...
#	@-grep --line-number --color 'assert(' $(mycontextdir)/*pp
#	@echo Find usage of 'cout' ...
	@-grep --line-number --color 'cout' $(mycontextdir)/*pp
#	@echo Find floating-point numbers ...
#	@-grep --line-number --color -E '\.*[0-9]' $(mycontextdir)/*pp
#	@-grep --line-number --color -E '(0-9)e' $(mycontextdir)/*pp
	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' $(mycontextdir)/*pp


########################################################################
# Target 'check' is a generic test. Currently, it defaults to 'tidy'

check: tidy

.PHONY: check


########################################################################
# Commands for cleaning out numerous files that are not part of the project.
# These remove the following:
# - clean: delete all binary files and temporary output, and the following two.
# - vtkclean: delete all *vtk output files
# - dependclean: delete .deps directories and their content

# TODO: rewrite the entire thing ...

# $(context).cleanpattern    := .all.o *.a *.o *.d *.so *.gch *.exe *.exe.stackdump *.out *.out.stackdump OUTPUT_CPPLINT.txt callgrind.out.* 
# $(context).vtkcleanpattern := *.vtk
# $(context).depcleanpattern := .deps

# $(context).cleanfiles    := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).cleanpattern   )) )
# $(context).vtkcleanfiles := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).vtkcleanpattern)) )
# $(context).depcleanfiles := $(patsubst %, $(CURDIR)/%, $(wildcard $($(context).depcleanpattern)) )

# cleanfiles    += $($(context).cleanfiles)
# vtkcleanfiles += $($(context).vtkcleanfiles)
# depcleanfiles += $($(context).depcleanfiles)

CMD_CLEAN    = rm -f .all.o *.a *.o *.d *.so *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_VTKCLEAN = rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk
CMD_DEPCLEAN = if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean vktclean dependclean
.PHONY: $(context).clean $(context).vktclean $(context).dependclean

clean:       $(context).clean
vtkclean:    $(context).vtkclean
dependclean: $(context).dependclean

$(context).clean $(context).vtkclean $(context).dependclean: mycontext    := $(context)
$(context).clean $(context).vtkclean $(context).dependclean: mycontextdir := $(contextdir)

$(context).clean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_CLEAN); $(CMD_VTKCLEAN); $(CMD_DEPCLEAN); 

$(context).vtkclean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_VTKCLEAN); 

$(context).dependclean: 
#	@-echo $(PWD)
	@-cd $(mycontextdir); $(CMD_DEPCLEAN); 


