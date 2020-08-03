
# Global definition of compiler and its flags 



# Do you want to use GCC or Clang?
# Comment out the appropriate definition below
FLAG_CXX := CLANG
# FLAG_CXX := GCC

# Do you want to ENABLE the use of tcmalloc?
# Comment out the following line to disable tcmalloc
# FLAG_USE_TCMALLOC=yes

# Do you want to ENABLE excessive warning options?
# Comment out the following line to disable excessive warning options
# FLAG_EXCESSIVE_WARNINGS=yes

# Do you want to ENABLE the use of openMP?
# Comment out the following line to disable compilation with openMP
# OPENMP_FLAG := -fopenmp

# Do you want to DISABLE checking of meshes?
# Comment out the following line to retain meaningful check routines for meshes
FLAG_DO_NOT_CHECK_MESHES := -DDO_NOT_CHECK_MESHES

# Do you want to enable static analysis during the compilation process
# Comment out the following line to disable static analysis
# FLAG_DO_STATICANALYSIS=yes

# Do you want to ENABLE the standard library debugging flags 
# Comment out the following line to disable the standard library debugging flags 
FLAG_DO_COMPILEDEBUGMODE := -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

# Do you want to DISABLE the general assert macro?
# Comment out the following line to disable the general assert macro
FLAG_DNDEBUG := -DNDEBUG

# Do you want to ENABLE the Clang sanitizer?
# Comment out the following line to disable compilation with the Clang sanitizer
# FLAG_DO_USE_SANITIZER=yes





### Compiler 

ifeq ($(FLAG_CXX),GCC)

  CXX := g++
  
else ifeq ($(FLAG_CXX),CLANG)

  CXX := clang++

else

  $(error No compiler recognized)

endif








### Language and dialects 

CXXFLAGS_LANG := -std=c++17 -pedantic 








### Format for diagnostic messages 

ifeq ($(FLAG_CXX),CLANG)

  CXXFLAGS_DIAGFORMAT := -fdiagnostics-show-template-tree

endif




### 
### TODO:
### Enable Wall extra pedantic 
### but filter out annoying warnings for gcc / clang 
### 
### Than add shared options and compiler specific options 
### 
### That can be found in the GCC man pages 
### 
### 


### Warnings 

CXXFLAGS_WARNINGS := 
CXXFLAGS_WARNINGS += 
CXXFLAGS_WARNINGS += -Wall -Wextra -Wpedantic -Wno-vla 
CXXFLAGS_WARNINGS += -Wno-conversion -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter

ifeq ($(FLAG_EXCESSIVE_WARNINGS),yes)

CXXFLAGS_WARNINGS += -Wodr -Wmissing-field-initializers -Wctor-dtor-privacy -Wsign-promo -Woverloaded-virtual -Wno-missing-braces

CXXFLAGS_WARNINGS += -Wundef 
CXXFLAGS_WARNINGS += -Wcast-align -Wmissing-declarations -Wredundant-decls -Wno-redundant-decls
CXXFLAGS_WARNINGS += -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wmisleading-indentation
# CXXFLAGS_WARNINGS += -Wold-style-cast

ifeq ($(FLAG_CXX),GCC)

  CXXFLAGS_WARNINGS +=  -Wmultiple-inheritance -Wvirtual-inheritance
  CXXFLAGS_WARNINGS += -Wpointer-arith
  CXXFLAGS_WARNINGS += -Wreturn-local-addr
  CXXFLAGS_WARNINGS += -Wcast-qual -Wcast-align
  CXXFLAGS_WARNINGS += -Wfloat-equal -Wdouble-promotion -Wfloat-conversion -Wno-type-limits
  CXXFLAGS_WARNINGS += -Wno-sign-compare
#   CXXFLAGS_WARNINGS += -Wconversion 
  CXXFLAGS_WARNINGS += -Wno-sign-conversion
  CXXFLAGS_WARNINGS += -Wlogical-op -Wreturn-local-addr
  #check which options there are ... 
  CXXFLAGS_WARNINGS += -Wnull-dereference
  CXXFLAGS_WARNINGS += -Wswitch-default -Wswitch-enum
  CXXFLAGS_WARNINGS += -Wunused -Wno-unused-variable
  CXXFLAGS_WARNINGS += -Wuninitialized
  CXXFLAGS_WARNINGS += -Wsuggest-final-types -Wsuggest-final-methods -Wsuggest-override
  CXXFLAGS_WARNINGS += -Walloca
  CXXFLAGS_WARNINGS += -Warray-bounds=2
  CXXFLAGS_WARNINGS += -Wduplicated-branches -Wduplicated-cond
  #CXXFLAGS_WARNINGS += -Wshadow=global
  CXXFLAGS_WARNINGS += -Wno-shadow
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += -Wdangling-else -Wparentheses
  CXXFLAGS_WARNINGS += -Wwrite-strings
  CXXFLAGS_WARNINGS += -Wunsafe-loop-optimizations -Winline -Wvector-operation-performance -Wdisabled-optimization
	#TODO What is stack smashing???? HSA????
	#TODO Read the format warnings 

else ifeq ($(FLAG_CXX),CLANG)

  CXXFLAGS_WARNINGS += -Wno-vla-extension -Werror-implicit -Wabsolute-value -Wno-shorten-64-to-32 -Walloca -Wanon-enum-enum-conversion

endif
 
endif
 




 
 
 
 
 
 ### Static analyzer 

ifeq ($(FLAG_DO_STATICANALYSIS),yes)

ifeq ($(FLAG_CXX),GCC)

CXXFLAGS_STATICANALYSER := -fanalyzer -Wanalyzer-too-complex

else ifeq ($(FLAG_CXX),CLANG)

CXXFLAGS_STATICANALYSER := 

endif

endif









### Debugging 

CXXFLAGS_DEBUG := -g








### Sanitizer instrumentation 

SANITIZER_FLAG := -pg -fno-omit-frame-pointer

ifeq ($(FLAG_DO_USE_SANITIZER),yes)

	SANITIZER_UNDEFINED_FLAG:=undefined,float-divide-by-zero,unsigned-integer-overflow,implicit-conversion,nullability-arg,nullability-assign,nullability-return

	# SANITIZER_SAFESTACK_FLAG:=safe-stack

	SANITIZER_ADDRESSLEAK_FLAG:=address,leak

	# SANITIZER_CFI_FLAG:=cfi

	# -fvisibility=hidden

	# Comment out the following line to disable ALL built-in sanitizers 
	SANITIZER_FLAG += -fsanitize=$(SANITIZER_UNDEFINED_FLAG),$(SANITIZER_ADDRESSLEAK_FLAG),$(SANITIZER_SAFESTACK_FLAG) -pg -fno-omit-frame-pointer
	
	SANITIZER_FLAG += -ftrapv
endif

ifeq ($(FLAG_USE_TCMALLOC),yes)
	CXXFLAGS_MALLOC=-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
else
	CXXFLAGS_MALLOC=
endif








# Optimization 

CXXFLAGS_OPTIMIZE := -Ofast
CXXFLAGS_OPTIMIZE += -march=native $(OPENMP_FLAG)
# CXXFLAGS_OPTIMIZE := -inline-threshold=1200
CXXFLAGS_OPTIMIZE += -flto








# Code generation options

CXXFLAGS_CODEGEN := -fpic 
CXXFLAGS_CODEGEN += -fno-exceptions -fvisibility=default










###############################################
#                                             #
#     CXXFLAGS - gather compiler flags        #
#                                             #
###############################################

CXXFLAGS := ${CXXFLAGS_LANG} ${CXXFLAGS_DIAGFORMAT} ${CXXFLAGS_WARNINGS} ${CXXFLAGS_STATICANALYSER} ${CXXFLAGS_DEBUG} $(SANITIZER_FLAG) ${CXXFLAGS_MALLOC} ${CXXFLAGS_OPTIMIZE} ${CXXFLAGS_CODEGEN}





###############################################
#                                             #
#            preprocessor flags               #
#                                             #
###############################################

CPPFLAGS := $(FLAG_DO_NOT_CHECK_MESHES) $(FLAG_DNDEBUG) $(FLAG_DO_COMPILEDEBUGMODE)



###############################################
#                                             #
#               linker flags                  #
#                                             #
###############################################

ifeq ($(FLAG_USE_TCMALLOC),yes)
	LDLIBS :=-l:libtcmalloc.so.4
else
	LDLIBS :=
endif




