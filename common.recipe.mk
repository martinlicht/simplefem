
# Global definition of compiler and its flags 



# Do you want to use GCC or Clang?
# Comment out the appropriate definition below
FLAG_CXX := CLANG
# FLAG_CXX := GCC

# Do you want to ENABLE the use of tcmalloc?
# Comment out the following line to disable tcmalloc
# FLAG_USE_TCMALLOC=yes

# Do you want to ENABLE the use of openMP?
# Comment out the following line to disable compilation with openMP
# OPENMP_FLAG := -fopenmp

# Do you want to DISABLE checking of meshes?
# Comment out the following line to retain meaningful check routines for meshes
# FLAG_DO_NOT_CHECK_MESHES := -DDO_NOT_CHECK_MESHES



# Do you want to ENABLE the debugging flags 
# Comment out the following line to disable the general debugging flags 
FLAG_DO_COMPILEDEBUGMODE := -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

# Do you want to DISABLE the general assert macro?
# Comment out the following line to disable the general assert macro
# FLAG_DNDEBUG := -DNDEBUG

# Do you want to ENABLE the Clang sanitizer?
# Comment out the following line to disable compilation with the Clang sanitizer
# FLAG_DO_USE_SANITIZER=yes





### Compiler 

ifeq ($(FLAG_CXX),GCC)

  CXX := g++
  
  CXXFLAGS_WARNINGS_COMPILER := -Wlogical-op -Wreturn-local-addr -Wunsafe-loop-optimizations -Wvector-operation-performance -Wno-type-limits

else ifeq ($(FLAG_CXX),CLANG)

  CXX := clang++

  CXXFLAGS_WARNINGS_COMPILER := -Werror-implicit -Wabsolute-value

else

  $(error No compiler recognized)

endif




### Language Specification 

CXXFLAGS_LANG := -std=c++17 -pedantic 





### Optimization flags 

CXXFLAGS_OPTIMIZE := -Ofast -march=native $(OPENMP_FLAG) 
#-mllvm -inline-threshold=1200
# -fopenmp 
# -O3 -Ofast -march=native -flto -frename-registers -frename-registers




### Compilation parameters 

ifeq ($(FLAG_DO_USE_SANITIZER),yes)
	SANITIZER_UNDEFINED_FLAG:=undefined,float-divide-by-zero,unsigned-integer-overflow,implicit-conversion,nullability-arg,nullability-assign,nullability-return
	# SANITIZER_SAFESTACK_FLAG:=safe-stack
	SANITIZER_ADDRESSLEAK_FLAG:=address,leak
	# SANITIZER_CFI_FLAG:=cfi
	# -fvisibility=hidden

	# Comment out the following line to disable ALL built-in sanitizers 
	SANITIZER_FLAG := -fsanitize=$(SANITIZER_UNDEFINED_FLAG),$(SANITIZER_ADDRESSLEAK_FLAG),$(SANITIZER_SAFESTACK_FLAG) -pg -fno-omit-frame-pointer
else
	SANITIZER_FLAG := -fsanitize=$(SANITIZER_UNDEFINED_FLAG),$(SANITIZER_ADDRESSLEAK_FLAG),$(SANITIZER_SAFESTACK_FLAG) -pg -fno-omit-frame-pointer
endif


CXXFLAGS_WARNINGS_EXTRA := -Wundef -Wcast-align -Wwrite-strings -Wmissing-declarations -Wredundant-decls -Woverloaded-virtual -Wformat -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wmisleading-indentation -Wunused -Wpedantic -Wpointer-arith -Wconversion -Wparentheses -Wdisabled-optimization -Wno-sign-conversion -Wno-shorten-64-to-32 -Wno-shadow

CXXFLAGS_WARNINGS := -Wall -Wextra -Wodr -Wno-unused-variable -Wno-sign-compare -Wno-missing-braces -Wmissing-field-initializers -Wfloat-equal $(CXXFLAGS_WARNINGS_EXTRA) $(CXXFLAGS_WARNINGS_COMPILER)


ifeq ($(FLAG_USE_TCMALLOC),yes)
	CXXFLAGS_MALLOC=-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
else
	CXXFLAGS_MALLOC=
endif

CXXFLAGS_MISC := -g -fpic -fno-exceptions

CXXFLAGS := ${CXXFLAGS_LANG} ${CXXFLAGS_WARNINGS} ${CXXFLAGS_OPTIMIZE} ${CXXFLAGS_MISC} $(SANITIZER_FLAG) ${CXXFLAGS_MALLOC}


### Preprocessor flags 

CPPFLAGS := $(FLAG_DO_NOT_CHECK_MESHES) $(FLAG_DNDEBUG) $(FLAG_DO_COMPILEDEBUGMODE)



### Linker Flags

ifeq ($(FLAG_USE_TCMALLOC),yes)
	LDLIBS :=-l:libtcmalloc.so.4
else
	LDLIBS :=
endif




