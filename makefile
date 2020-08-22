SHELL = /bin/sh

default: all

all: lib test 

buildobjects:
	cd ./utility && $(MAKE) buildobjects
	cd ./combinatorics && $(MAKE) buildobjects 
	cd ./operators && $(MAKE) buildobjects
	cd ./dense && $(MAKE) buildobjects
	cd ./sparse && $(MAKE) buildobjects
	cd ./solver && $(MAKE) buildobjects
	cd ./mesh && $(MAKE) buildobjects
	cd ./vtk && $(MAKE) buildobjects
	cd ./matrixmarket && $(MAKE) buildobjects
	# cd ./fem && $(MAKE)
	echo "finished building objects" 

buildso:
	cd ./utility && $(MAKE) buildso
	cd ./combinatorics && $(MAKE) buildso
	cd ./operators && $(MAKE) buildso
	cd ./dense && $(MAKE) buildso
	cd ./sparse && $(MAKE) buildso
	cd ./solver && $(MAKE) buildso
	cd ./mesh && $(MAKE) buildso
	cd ./vtk && $(MAKE) buildso
	cd ./matrixmarket && $(MAKE) buildso
	# cd ./fem && $(MAKE)
	echo "finished building shared libraries." 

lib:
	cd ./utility && $(MAKE) 
	cd ./combinatorics && $(MAKE) 
	cd ./operators && $(MAKE) 
	cd ./dense && $(MAKE) 
	cd ./sparse && $(MAKE) 
	cd ./solver && $(MAKE) 
	cd ./mesh && $(MAKE) 
	cd ./vtk && $(MAKE)
	cd ./matrixmarket && $(MAKE)
	# cd ./fem && $(MAKE)
	echo "finished building objects and shared libraries." 


.PHONY: test 
test:
	cd ./tests && $(MAKE)
	echo "finished building tests."


.PHONY: vtkclean
vtkclean:
	-rm -f *.vtk
	-rm -f ./*/*.vtk
	-rm -f ./*/*/*.vtk


.PHONY: clean
clean: 
	cd ./utility && $(MAKE) clean
	cd ./combinatorics && $(MAKE) clean
	cd ./operators && $(MAKE) clean
	cd ./dense && $(MAKE) clean
	cd ./sparse && $(MAKE) clean
	cd ./solver && $(MAKE) clean
	cd ./mesh && $(MAKE) clean
	cd ./vtk && $(MAKE) clean
	cd ./matrixmarket && $(MAKE) clean
	cd ./fem && $(MAKE) clean
	cd ./tests && $(MAKE) clean
	$(MAKE) -f common.upkeep.mk clean
	@echo "finished cleaning." 


.PHONY: dependclean
dependclean: 
	cd ./utility && $(MAKE) dependclean
	cd ./combinatorics && $(MAKE) dependclean
	cd ./operators && $(MAKE) dependclean
	cd ./dense && $(MAKE) dependclean
	cd ./sparse && $(MAKE) dependclean
	cd ./solver && $(MAKE) dependclean
	cd ./mesh && $(MAKE) dependclean
	cd ./vtk && $(MAKE) dependclean
	cd ./matrixmarket && $(MAKE) dependclean
	cd ./fem && $(MAKE) dependclean
	cd ./tests && $(MAKE) dependclean
	$(MAKE) -f common.upkeep.mk dependclean
	@echo "finished cleaning dependencies." 


.PHONY: tidy
tidy: 
	cd ./utility && $(MAKE) tidy
	cd ./combinatorics && $(MAKE) tidy
	cd ./operators && $(MAKE) tidy
	cd ./dense && $(MAKE) tidy
	cd ./sparse && $(MAKE) tidy
	cd ./solver && $(MAKE) tidy
	cd ./mesh && $(MAKE) tidy
	cd ./vtk && $(MAKE) tidy
	cd ./matrixmarket && $(MAKE) tidy
	cd ./fem && $(MAKE) tidy
	cd ./tests && $(MAKE) tidy
	$(MAKE) -f makefile.tidy tidy
	@echo "finished tidying." 


.PHONY: check
check: 
	cd ./utility && $(MAKE) check
	cd ./combinatorics && $(MAKE) check
	cd ./operators && $(MAKE) check
	cd ./dense && $(MAKE) check
	cd ./sparse && $(MAKE) check
	cd ./solver && $(MAKE) check
	cd ./mesh && $(MAKE) check
	cd ./vtk && $(MAKE) check
	cd ./matrixmarket && $(MAKE) check
	cd ./fem && $(MAKE) check
# 	cd ./tests && $(MAKE) check
	$(MAKE) -f makefile.check check
	@echo "finished checking." 


# shared: default
# 	mkdir build
# 	gcc -shared -o ./build/combinatorics.so ./combinatorics/*.o
# 	gcc -shared -o ./build/operators.so ./operators/*.o
# 	gcc -shared -o ./build/dense.so ./dense/*.o
# 	gcc -shared -o ./build/sparse.so ./sparse/*.o
# 	gcc -shared -o ./build/solver.so ./solver/*.o
# 	gcc -shared -o ./build/mesh.so ./mesh/*.o
# 	gcc -shared -o ./build/vtk.so ./vtk/*.o
# 
# cleanshared:
# 	rm -f ./build/*.so




CHECK_OPTION= --enable=warning,style,performance,portability --std=c++11 -q
CHECK_FILES= basic.hpp basic/*.?pp combinatorics/*.?pp operators/*.?pp


# .PHONY: check
# check:
# 	cppcheck $(CHECK_OPTION) . -ilegacy/ -iplayground/
# 	cppcheck $(CHECK_OPTION) basic.hpp
# 	cppcheck $(CHECK_OPTION) basic/*.?pp
# 	cppcheck $(CHECK_OPTION) combinatorics/*.?pp
# 	cppcheck $(CHECK_OPTION) operators/*.?pp
# 	cppcheck $(CHECK_OPTION) dense/*.?pp
# 	cppcheck $(CHECK_OPTION) sparse/*.?pp
# 	cppcheck $(CHECK_OPTION) solver/*.?pp
# 	cppcheck $(CHECK_OPTION) mesh/*.?pp
# 	cppcheck $(CHECK_OPTION) */*.?pp



# 	cd ./utility && $(MAKE) cpplint
# 	cd ./combinatorics && $(MAKE) cpplint
# 	cd ./operators && $(MAKE) cpplint
# 	cd ./dense && $(MAKE) cpplint
# 	cd ./sparse && $(MAKE) cpplint
# 	cd ./solver && $(MAKE) cpplint
# 	cd ./mesh && $(MAKE) cpplint
# 	cd ./vtk && $(MAKE) cpplint
# 	cd ./matrixmarket && $(MAKE) cpplint
# 	cd ./fem && $(MAKE) cpplint
# # 	cd ./tests && $(MAKE) cpplint
# 	$(MAKE) -f makefile.cpplint cpplint
# 	@echo "finished cpplinting." 


.PHONY: cpplint
cpplint:
	( ./cpplint.py --exclude=tests/* --exclude=tests/*/* --exclude=.legacy/* --exclude=.private/* --exclude=.playground/* --recursive --filter=-whitespace,-legal --quiet . ) | sort | uniq -c > OUTPUT_CPPLINT.txt



