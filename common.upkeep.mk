
# This file contains a collection of common commands 
# for the 'upkeep' of the project directories.
# In particular, cleaning off temporary files and static code analysis,
# which can be executed in each single source directory


# clean out executables, debug output, and intermediate files

.PHONY: clean
clean: vtkclean dependclean
# 	@echo clean content ...
#	@-rm -f OUTPUT_CPPLINT.txt
#	@-rm -f callgrind.out.*
#	@-rm -f .all.o *.a *.o *.d *.so *.gch
#	@-rm -f *.exe *.exe.stackdump
#	@-rm -f *.out *.out.stackdump 
	@-rm -f \
	 .all.o *.a *.o *.d *.so *.gch \
	 OUTPUT_CPPLINT.txt \
	 callgrind.out.* \
	 *.exe *.exe.stackdump \
	 *.out *.out.stackdump 


# clean out and remove the dependency directories 

.PHONY: dependclean
dependclean:
# 	@echo clean dependency files
	@-if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; rmdir .deps/; fi 
#	@-if [ -d .deps/ ]; then rmdir .deps/; fi 



# clean out all VTK files

.PHONY: vtkclean
vtkclean:
	@-rm -f ./*.vtk ./*/*.vtk ./*/*/*.vtk


# apply clang-tidy to all cpp and hpp files in the directory

.PHONY: tidy
tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++2a


# apply cppcheck to all cpp and hpp files in the directory 

.PHONY: cppcheck
cppcheck:
	cppcheck -i ./.playground/ -i ./.legacy --enable=warning,style,performance,portability --suppress=duplicateCondition --suppress=assertWithSideEffect --suppress=useStlAlgorithm --std=c++17 -q . ./*pp


# apply cpplint to all cpp and hpp files in the directory 

.PHONY: cpplint
cpplint:
	$(error This command is not implemented.)
	( ./../Tools/cpplint.py --exclude=.private/ --exclude=.legacy/ --exclude=.playground/ --recursive --filter=-whitespace,-legal,-build/namespace,readability/alt_tokens,readability/todo,readability/inheritance --quiet . ) | sort | uniq -c 2> OUTPUT_CPPLINT.txt


# regex several useful things 
# - find trailing white spaces 
# - find non-ASCII characters 
# - find consecutive spaces 

.PHONY: grepissues
grepissues:
	@echo Find trailing whitespace...
#	@-grep --line-number --color '\s+$$' -r ./*pp
#	@echo Find non-ASCII characters...
#	@-grep --line-number --color '[^\x00-\x7F]' -r ./*pp
#	@echo Find consecutive spaces...
#	@-grep --line-number --color '\b\s{2,}' -r ./*pp
#	@-grep --line-number --color 'assert(' ./*pp
	@-grep --line-number --color 'cout' ./*pp
#	@-grep --line-number --color -E '\.*[0-9]' ./*pp
	@-grep --line-number --color -E '(0-9)e' ./*pp
	@-grep --line-number --color -E '([0-9]+e[0-9]+)|([0-9]+\.[0-9]+)|((+-\ )\.[0-9]+)|((+-\ )[0-9]+\.)' ./*pp

