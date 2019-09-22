
# Global definition of compiler and its flags 

### Compiler 

CXX = clang++-8
# CXX = g++


### Compilation parameters 

CXXFLAGS_LANG = -std=c++2a -pedantic 

CXXFLAGS_WARNINGS = -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit

CXXFLAGS_OPTIMIZE = -Ofast -march=native -fopenmp
# -fopenmp 
# -O3 -Ofast -march=native -flto -frename-registers -frename-registers

CXXFLAGS_MISC =  -g -fpic -fno-exceptions -fsanitize=leak,address,undefined 
#-fsanitize=address,leak,undefined -pg -fno-omit-frame-pointer

CXXFLAGS = ${CXXFLAGS_LANG} ${CXXFLAGS_WARNINGS} ${CXXFLAGS_OPTIMIZE} ${CXXFLAGS_MISC}


### Preprocessor flags 

CPPFLAGS = 
#-DNDEBUG -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC




# *.o: ../common.make ./makefile




# .PHONY: tidy
# tidy:
# 	clang-tidy ./*.?pp -checks=llvm*,bugprone-*,clang-analyzer-*,misc-*,-llvm-header-guard,-llvm-include-order -- -std=c++17

# .PHONY: check
# check:
# 	cppcheck --enable=warning,style,performance,portability --std=c++17 -q .

# .PHONY: cpplint
# cpplint:
# 	( ./../cpplint.py --recursive --filter=-whitespace,-legal --quiet . ) | sort | uniq -c > TESTFOO





# CXXFLAGS = -std=c++17 -O3 -g -fno-exceptions -fpic -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit

# https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
 

# CXX = g++ -O0 -g -std=c++11 -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit
# CXX = clang++-6.0 -O0 -g -std=c++17 -fno-exceptions -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit
# CXX =     g++     -O0 -g -fno-exceptions -std=c++1z -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit
# CXX =     g++     -O0 -g -fsanitize=leak -std=c++1z -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit
# CXX =     g++ -O0 -g -fno-exceptions -std=c++14 -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -pedantic -Wall -Wextra -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Werror=implicit



