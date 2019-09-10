SHELL = /bin/sh

default:
	cd ./combinatorics && $(MAKE) 
	cd ./operators && $(MAKE) 
	cd ./dense && $(MAKE) 
	cd ./sparse && $(MAKE) 
	cd ./solver && $(MAKE) 
	cd ./mesh && $(MAKE) 
	cd ./vtk && $(MAKE)
	cd ./matrixmarket && $(MAKE)
	cd ./fem && $(MAKE)
	cd ./tests && $(MAKE)
	echo "finished building." 

clean: 
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
	$(MAKE) -f makefile.clean clean
	@echo "finished cleaning." 

tidy: 
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

check:
	cppcheck $(CHECK_OPTION) . -ilegacy/ -iplayground/
# 	cppcheck $(CHECK_OPTION) basic.hpp
# 	cppcheck $(CHECK_OPTION) basic/*.?pp
# 	cppcheck $(CHECK_OPTION) combinatorics/*.?pp
# 	cppcheck $(CHECK_OPTION) operators/*.?pp
# 	cppcheck $(CHECK_OPTION) dense/*.?pp
# 	cppcheck $(CHECK_OPTION) sparse/*.?pp
# 	cppcheck $(CHECK_OPTION) solver/*.?pp
# 	cppcheck $(CHECK_OPTION) mesh/*.?pp
# 	cppcheck $(CHECK_OPTION) */*.?pp

cpplint:
	( ./cpplint.py --recursive --filter=-whitespace,-legal --quiet . --exclude=tests/* --exclude=tests/*/* --exclude=legacy/* --exclude=playground/* ) | sort | uniq -c > TESTFOO



