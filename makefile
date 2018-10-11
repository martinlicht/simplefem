

default:
	cd ./combinatorics && $(MAKE) 
	cd ./operators && $(MAKE) 
	cd ./dense && $(MAKE) 
	cd ./sparse && $(MAKE) 
	cd ./solver && $(MAKE) 
	cd ./mesh && $(MAKE) 
	cd ./vtk && $(MAKE)
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
	cd ./tests && $(MAKE) clean
	$(MAKE) -f $(MAKE)file.clean clean
	@echo "finished cleaning." 

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



