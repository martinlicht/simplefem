
# This file contains a collection of common commands 
# for the 'upkeep' of the project directories.
# In particular, cleaning off temporary files and static code analysis,
# which can be executed in each single source directory



.PHONY: tidy
tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++17

.PHONY: check
check:
	cppcheck --enable=warning,style,performance,portability --suppress=duplicateCondition --suppress=assertWithSideEffect --suppress=useStlAlgorithm --std=c++17 -q . ./*/*pp

.PHONY: cpplint
cpplint:
	( ./../cpplint.py --exclude=.private/ --exclude=.legacy/ --exclude=.playground/ --recursive --filter=-whitespace,-legal,-build/namespace --quiet . ) | sort | uniq -c > OUTPUT_CPPLINT.txt

.PHONY: vtkclean
vtkclean:
	-rm -f *.vtk
	-rm -f ./*/*.vtk
	-rm -f ./*/*/*.vtk

.PHONY: dependclean
dependclean:
	-if [ -d .deps/ ]; then rm -f .deps/*.d .deps/.all.d; fi 
	-if [ -d .deps/ ]; then rmdir .deps/; fi 

.PHONY: clean
clean: vtkclean dependclean
	-rm -f OUTPUT_CPPLINT.txt
	-rm -f callgrind.out.*
	-rm -f .all.o *.o *.d *.so *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 

