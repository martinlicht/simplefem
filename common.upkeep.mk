
# This file contains a collection of common commands 
# for the 'upkeep' of the project directories.
# In particular, cleaning off temporary files and static code analysis,
# which can be executed in each single source directory


# clean out executables, debug output, and intermediate files

.PHONY: clean
clean: vtkclean dependclean
# 	@echo clean content ...
	@-rm -f OUTPUT_CPPLINT.txt
	@-rm -f callgrind.out.*
	@-rm -f .all.o *.a *.o *.d *.so *.gch
	@-rm -f *.exe *.exe.stackdump
	@-rm -f *.out *.out.stackdump 


# clean out and remove the dependency directories 

.PHONY: dependclean
dependclean:
# 	@echo clean dependency files
	@-if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; fi 
	@-if [ -d .deps/ ]; then rmdir .deps/; fi 



# clean out all VTK files

.PHONY: vtkclean
vtkclean:
	@-rm -f ./*.vtk
	@-rm -f ./*/*.vtk
	@-rm -f ./*/*/*.vtk


# apply clang-tidy to all cpp and hpp files in the directory

.PHONY: tidy
tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++17


# apply cppcheck to all cpp and hpp files in the directory 

.PHONY: cppcheck
cppcheck:
	cppcheck --enable=warning,style,performance,portability --suppress=duplicateCondition --suppress=assertWithSideEffect --suppress=useStlAlgorithm --std=c++17 -q . ./*pp


# apply cpplint to all cpp and hpp files in the directory 

.PHONY: cpplint
cpplint:
	( ./../Tools/cpplint.py --exclude=.private/ --exclude=.legacy/ --exclude=.playground/ --recursive --filter=-whitespace,-legal,-build/namespace,readability/alt_tokens,readability/todo,readability/inheritance --quiet . ) | sort | uniq -c > OUTPUT_CPPLINT.txt

