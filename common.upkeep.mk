
.PHONY: tidy
tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++17

.PHONY: check
check:
	cppcheck --enable=warning,style,performance,portability --suppress=duplicateCondition --suppress=assertWithSideEffect --suppress=useStlAlgorithm --std=c++17 -q . ./*/*pp

# .PHONY: cpplint
# cpplint:
# 	( ./../cpplint.py --exclude=.private/ --exclude=.legacy/ --exclude=.playground/ --recursive --filter=-whitespace,-legal,-build/namespace --quiet . ) | sort | uniq -c > OUTPUT_CPPLINT.txt

.PHONY: vtkclean
vtkclean:
	-rm -f *.vtk
	-rm -f ./*/*.vtk
	-rm -f ./*/*/*.vtk

.PHONY: clean
clean: vtkclean
	-if [ -d .deps/ ]; then rm -f .deps/*.d; fi 
	-if [ -d .deps/ ]; then rmdir .deps/; fi 
	-rm -f *.o *.d *.so *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 

.PHONY: dependclean
dependclean:
	-if [ -d .deps/ ]; then rm -f .deps/*.d; fi 
	-if [ -d .deps/ ]; then rmdir .deps/; fi 

