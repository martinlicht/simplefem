

default:
	cd ./combinatorics && make && cd ..
	cd ./operators && make && cd ..
	cd ./dense && make && cd ..
	cd ./sparse && make && cd ..
	cd ./solver && make && cd ..
	cd ./mesh && make && cd ..
	cd ./vtk && make && cd ..
	cd ./tests && make && cd ..
	echo "finished building." 

clean:
	cd ./combinatorics && make clean && cd ..
	cd ./operators && make clean && cd ..
	cd ./dense && make clean && cd ..
	cd ./sparse && make clean && cd ..
	cd ./solver && make clean && cd ..
	cd ./mesh && make clean && cd ..
	cd ./vtk && make clean && cd ..
	cd ./tests && make clean && cd ..
	make -f makefile.clean clean
	echo "finished cleaning." 

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







