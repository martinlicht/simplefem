
# This file contains makefile rules for object creation
# that can be applied in every single source directory
# They create the dependency auxiliary files, 
# the object files, and the shared libraries. 


# Do we need basename or path here?
# https://www.gnu.org/software/make/manual/html_node/File-Name-Functions.html

$(module).depdir  := $(moddir)/.deps

$(module).sources      := $(wildcard $(moddir)/*.cpp)
$(module).headers      := $(wildcard $(moddir)/*.hpp)
$(module).objects      := $(patsubst %.cpp,%.o,$($(module).sources))
$(module).dependencies := $(patsubst %.cpp,$($(module).depdir)/%.d,$(notdir $($(module).sources)))

$(module).sharedlibrarybasename := $(moddir)/lib$(module)
$(module).libraryobject         := $(moddir)/lib$(module).o
$(module).sharedlibrary         := $(moddir)/lib$(module).so
$(module).staticlibrary         := $(moddir)/lib$(module).a

# TODO Remove that idea with build here, and reduce the number of variables expected here. 
#      Instead, just produce $(module).build, and let the outside take care of that target. 

# TODO Can we introduce the recipe variables for all $(module).* components simultaneously?
$(module).%: mymodule := $(module)
$(module).%: mymoddir := $(moddir)

###################################################################################################
# All object files that are compiled (and their descendants) also depend on the makefiles 

$(moddir)/*.o $(moddir)/.all.o: $(moddir)/makefile
$(moddir)/*.o $(moddir)/.all.o: $(projectdir)/makefile $(projectdir)/common.compile.mk $(projectdir)/common.module.mk 


###################################################################################################
# How to compile the .o files and .all.o file 


$($(module).depdir):
	@echo $@
	@"mkdir" -p $@

DEPFLAGS = -MT $@ -MMD -MP -MF $($(mymodule).depdir)/$(notdir $*.d)
# Generate the dependency files 
# -MT $@ : sets the target of the makefile rule, stripped of any directory components
# -MP    : add dependencies as phony targets
# -MF ...: sets the output file for the rules 
# -MMD   : list headers as a by-product of compiling, excluding system headers
# alternatively,
# -MM    : list headers, excluding system headers, and do not compile anything


.PHONY: make_dependencies $(module).make_dependencies
make_dependencies: $(module).make_dependencies
# $(module).make_dependencies: mymodule := $(module)
$(module).make_dependencies: $($(module).depdir)
	@for item in $($(mymodule).sources); do \
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $$item -MM -MP -MF $($(module).depdir)/$$item.d; \
	done


$(moddir)/$($(module).objects): %.o: %.cpp $($(module).depdir)/%.d | $($(module).depdir)
	@echo Compiling and listing dependencies: $@ 
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS)                       $< -c -o $@ $(DEPFLAGS)

# $(module).all.o: mymodule := $(module)
# $(module).all.o: mymoddir := $(moddir)
$(moddir)/.all.o: $($(module).sources) $(moddir)/.all.cpp $($(module).depdir)/.all.d | $($(module).depdir)
	@echo $(mymodule).depdir
	@echo $($(mymodule).depdir) 
	@echo $($(mymodule).sources) 
	@echo $($(mymodule).headers) 
	@echo $($(mymodule).objects) 
	@echo $($(mymodule).dependencies) 
	@echo $($(mymodule).sharedlibrarybasename)
	@echo $($(mymodule).libraryobject)
	@echo $($(mymodule).sharedlibrary)
	@echo $($(mymodule).staticlibrary)
	@echo $@
	@echo $*
	@echo Compiling and listing dependencies: $($(mymodule).libraryobject)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(mymoddir)/.all.cpp -c -o $@  $(DEPFLAGS)

$($(module).depdir)/.all.d:

$($(module).dependencies):

-include $($(module).depdir)/.all.d
-include $($(module).dependencies)


# How to provide the library .o file, and possibly static/shared libraries

# $($(module).libraryobject): mymodule := $(module)
# $($(module).libraryobject): mymoddir := $(moddir)
$($(module).libraryobject): $(moddir)/.all.o
	@echo Library object: $@
	@cp $(mymoddir)/.all.o $($(mymodule).libraryobject)

# $($(module).sharedlibrary): mymodule := $(module)
# $($(module).sharedlibrary): mymoddir := $(moddir)
$($(module).sharedlibrary): $($(module).libraryobject)
	@echo Shared library: $@
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(LDLIBS)

# $($(module).staticlibrary): mymodule := $(module)
# $($(module).staticlibrary): mymoddir := $(moddir)
$($(module).staticlibrary): $($(module).libraryobject)
	@echo Static library: $@
	@ar rcs $($(mymodule).staticlibrary) $($(mymodule).libraryobject)



###################################################################################################
# Finally, determine the build target depending on whether shared libraries are poossible or not

.PHONY: buildobjects buildso builda
.PHONY: $(module).buildobjects $(module).buildso $(module).builda

$(module).buildobjects: $($(module).objects)
$(module).buildso:      $($(module).sharedlibrary)
$(module).builda:       $($(module).staticlibrary)

buildobjects: $(module).buildobjects
buildso:      $(module).buildso
builda:       $(module).builda


.PHONY: $(buildtarget) $(module).build

$(buildtarget): $(module).build

ifeq ($(OS),Windows_NT)
$(module).build: $(module).builda
else
$(module).build: $(module).builda $(module).buildso
endif



########################################################################
# List several objects. Read-only

.PHONY:  list_of_objects $(module).list_of_objects
.SILENT: list_of_objects $(module).list_of_objects
list_of_objects: $(module).list_of_objects
# $(module).list_of_objects: mymodule := $(module)
# $(module).list_of_objects: mymoddir := $(moddir)
$(module).list_of_objects: 
	@echo $(projecdir);
	@echo $(mymodule);
	@echo $(mymoddir);
	@echo $($(mymodule).staticlibrary);
	@echo $($(mymodule).sharedlibrary);
	@echo $($(mymodule).objects);
	@echo $($(mymodule).dependencies);
	@echo $($(mymodule).build);


########################################################################
# Check whether the header files have correct syntax. Read-only.

$(module).headerchecks := $(patsubst %.hpp,check-%.hpp,$($(module).headers))

.PHONY: checkheaders $($(module).headerchecks)
checkheaders: $($(module).headerchecks)
$($(module).headerchecks): check-%.hpp : 
	$(info Check header: $*.hpp)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $*.hpp -fsyntax-only




########################################################################
# Apply clang-tidy to all cpp and hpp files in the directory. Read-only.

.PHONY: tidy $(module).tidy
tidy: $(module).tidy
$(module).tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++2a


########################################################################
# Apply cppcheck to all cpp and hpp files in the directory. Read-only. 

.PHONY: cppcheck $(module).cppcheck
cppcheck: $(module).cppcheck
$(module).cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy \
	--enable=warning,style,performance,portability --suppress=duplicateCondition\
	--suppress=assertWithSideEffect --suppress=useStlAlgorithm\
	--std=c++17 -q . ./*pp


########################################################################
# Regex several useful things. Read-only. 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues $(module).grepissues
grepissues: $(module).grepissues
# $(module).grepissues: mymodule := $(module)
# $(module).grepissues: mymoddir := $(moddir)
$(module).grepissues:
#	@echo Find trailing whitespace...
#	@-grep --line-number --color '\s+$$' -r $(mymoddir)/*pp
#	@echo Find non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r $(mymoddir)/*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r $(mymoddir)/*pp
#	@-grep --line-number --color 'assert(' $(mymoddir)/*pp
	@-grep --line-number --color 'cout' $(mymoddir)/*pp
#	@-grep --line-number --color -E '\.*[0-9]' $(mymoddir)/*pp
	@-grep --line-number --color -E '(0-9)e' $(mymoddir)/*pp
	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' $(mymoddir)/*pp


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

CMD_CLEAN    = rm -f .all.o *.a *.o *.d *.so *.gch OUTPUT_CPPLINT.txt callgrind.out.* *.exe *.exe.stackdump *.out *.out.stackdump 
CMD_VTKCLEAN = rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk
CMD_DEPCLEAN = if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 

.PHONY: clean vktclean dependclean
.PHONY: $(module).clean $(module).vktclean $(module).dependclean

clean:       $(module).clean
vtkclean:    $(module).vtkclean
dependclean: $(module).dependclean

$(module).clean $(module).vtkclean $(module).dependclean: mymodule := $(module)
$(module).clean $(module).vtkclean $(module).dependclean: mymoddir := $(moddir)

$(module).clean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_CLEAN); $(CMD_VTKCLEAN); $(CMD_DEPCLEAN); 

$(module).vtkclean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_VTKCLEAN); 

$(module).dependclean: 
#	@-echo $(PWD)
	@-cd $(mymoddir); $(CMD_DEPCLEAN); 


