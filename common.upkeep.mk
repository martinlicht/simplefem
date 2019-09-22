
.PHONY: tidy
tidy:
	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++17

.PHONY: check
check:
	cppcheck --enable=warning,style,performance,portability --std=c++17 -q .

.PHONY: cpplint
cpplint:
	( ./../cpplint.py --recursive --filter=-whitespace,-legal --quiet . ) | sort | uniq -c > TESTFOO


clean: 
	-rm -f .deps/*.d
	-rm -f *.o *.d *.so *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 
