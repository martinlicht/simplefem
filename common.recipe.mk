
# Global definition of compiler and its flags 



# Do you want to use GCC or Clang?
# Uncomment the appropriate definition below
FLAG_CXX := CLANG
# FLAG_CXX := GCC

# Do you want to ENABLE the use of tcmalloc?
# Uncomment the following line to enable tcmalloc
# FLAG_USE_TCMALLOC=yes

# Do you want to ENABLE excessive warning options?
# Uncomment the following line to enable excessive warning options
# FLAG_EXCESSIVE_WARNINGS=yes

# Do you want to ENABLE the use of openMP?
# Uncomment the following line to enable compilation with openMP
# OPENMP_FLAG := -fopenmp

# Do you want to DISABLE checking of meshes?
# Uncomment the following line to disable extensive check routines for meshes
FLAG_DO_NOT_CHECK_MESHES := -DDO_NOT_CHECK_MESHES

# Do you want to enable static analysis during the compilation process
# Uncomment the following line to enable static analysis
# FLAG_DO_STATICANALYSIS=yes

# Do you want to ENABLE the standard library debugging flags 
# Uncomment the following line to enable the standard library debugging flags 
FLAG_DO_COMPILEDEBUGMODE := -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC

# Do you want to DISABLE the general assert macro?
# Uncomment the following line to disable the general assert macro
# FLAG_DONTASSERT := -DNDEBUG

# Do you want to ENABLE the Clang sanitizer?
# Uncomment the following line to enable compilation with the Clang sanitizer
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
CXXFLAGS_WARNINGS += -Wall -Wextra -Wpedantic

ifeq ($(FLAG_EXCESSIVE_WARNINGS),yes)

CXXFLAGS_WARNINGS += -Wodr -Wmissing-field-initializers -Wctor-dtor-privacy -Wsign-promo -Woverloaded-virtual -Wno-missing-braces

CXXFLAGS_WARNINGS += -Wundef 
CXXFLAGS_WARNINGS += -Wcast-align -Wcast-qual -Wmissing-declarations -Wredundant-decls -Wno-redundant-decls
CXXFLAGS_WARNINGS += -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wmisleading-indentation
# CXXFLAGS_WARNINGS += -Wold-style-cast

ifeq ($(FLAG_CXX),GCC)

  CXXFLAGS_WARNINGS +=  -Wmultiple-inheritance -Wvirtual-inheritance
  CXXFLAGS_WARNINGS += -Wpointer-arith
  CXXFLAGS_WARNINGS += -Wreturn-local-addr
  CXXFLAGS_WARNINGS += -Wfloat-equal -Wdouble-promotion -Wfloat-conversion
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

  CXXFLAGS_WARNINGS += -Wabstract-vbase-init
  CXXFLAGS_WARNINGS += -Walloca
  CXXFLAGS_WARNINGS += -Wno-vla-extension -Werror-implicit -Wabsolute-value -Wno-shorten-64-to-32
  CXXFLAGS_WARNINGS += -Wanon-enum-enum-conversion -Wassign-enum
  CXXFLAGS_WARNINGS += -Warray-bounds-pointer-arithmetic
  CXXFLAGS_WARNINGS += -Wbad-function-cast
  CXXFLAGS_WARNINGS += -Wc++11-narrowing
  CXXFLAGS_WARNINGS += -Wchar-subscripts
  CXXFLAGS_WARNINGS += -Wclass-varargs
  CXXFLAGS_WARNINGS += -Wcomma
  CXXFLAGS_WARNINGS += -Wcomment
  CXXFLAGS_WARNINGS += -Wconsumed
  CXXFLAGS_WARNINGS += -Wconversion -Wno-sign-conversion
  CXXFLAGS_WARNINGS += -Wctad-maybe-unsupported
  CXXFLAGS_WARNINGS += -Wdate-time
  CXXFLAGS_WARNINGS += -Wdelete-non-abstract-non-virtual-dtor
  CXXFLAGS_WARNINGS += -Wdeprecated
  #CXXFLAGS_WARNINGS += -Wdouble-promotion
  #CXXFLAGS_WARNINGS += -Wdtor-name  TODO: Is this one actually defined???
  CXXFLAGS_WARNINGS += -Wduplicate-decl-specifier -Wduplicate-enum -Wduplicate-method-arg -Wduplicate-method-match 
  CXXFLAGS_WARNINGS += -Wembedded-directive
  CXXFLAGS_WARNINGS += -Wempty-init-stmt
  CXXFLAGS_WARNINGS += -Wenum-compare-conditional
  CXXFLAGS_WARNINGS += -Wexceptions
  CXXFLAGS_WARNINGS += -Wextra-semi
  CXXFLAGS_WARNINGS += -Wfloat-conversion -Wfloat-equal
  CXXFLAGS_WARNINGS += -Wformat
  CXXFLAGS_WARNINGS += -Wheader-hygiene
  CXXFLAGS_WARNINGS += -Widiomatic-parentheses
  CXXFLAGS_WARNINGS += -Wimplicit-float-conversion
  CXXFLAGS_WARNINGS += -Wimplicit-function-declaration
  CXXFLAGS_WARNINGS += -Winfinite-recursion
  CXXFLAGS_WARNINGS += -Wint-conversion
  CXXFLAGS_WARNINGS += -Wkeyword-macro
  #CXXFLAGS_WARNINGS += -Wmain
  CXXFLAGS_WARNINGS += -Wmisleading-indentation
  CXXFLAGS_WARNINGS += -Wmissing-braces
  CXXFLAGS_WARNINGS += -Wmissing-field-initializers
  CXXFLAGS_WARNINGS += -Wmissing-prototypes
  CXXFLAGS_WARNINGS += -Wmissing-variable-declarations
  CXXFLAGS_WARNINGS += -Wmost -Wmove 
  CXXFLAGS_WARNINGS += -Wnewline-eof
  CXXFLAGS_WARNINGS += -Wnon-virtual-dtor
  CXXFLAGS_WARNINGS += -Wnonportable-system-include-path
  CXXFLAGS_WARNINGS += -Wnull-pointer-arithmetic 
  CXXFLAGS_WARNINGS += -Woverlength-strings
  CXXFLAGS_WARNINGS += -Woverloaded-virtual
  CXXFLAGS_WARNINGS += -Winitializer-overrides
  CXXFLAGS_WARNINGS += -Woverriding-method-mismatch
  CXXFLAGS_WARNINGS += -Wparentheses
  CXXFLAGS_WARNINGS += -Wpessimizing-move
  CXXFLAGS_WARNINGS += -Wrange-loop-analysis
  CXXFLAGS_WARNINGS += -Wredundant-move
  CXXFLAGS_WARNINGS += -Wredundant-parens
#   CXXFLAGS_WARNINGS += -Wredundant-parentheses
  CXXFLAGS_WARNINGS += -Wreserved-id-macro
  CXXFLAGS_WARNINGS += -Wreserved-user-defined-literal
  CXXFLAGS_WARNINGS += -Wreturn-std-move
  CXXFLAGS_WARNINGS += -Wself-assign -Wself-move
  CXXFLAGS_WARNINGS += -Wsemicolon-before-method-body
  CXXFLAGS_WARNINGS += -Wsometimes-uninitialized
  CXXFLAGS_WARNINGS += -Wstrict-prototypes
  CXXFLAGS_WARNINGS += -Wstring-conversion
#   CXXFLAGS_WARNINGS += -Wsuggest-destructor-override
  CXXFLAGS_WARNINGS += -Wtautological-compare
  CXXFLAGS_WARNINGS += -Wtautological-type-limit-compare
#   CXXFLAGS_WARNINGS += -Wuninitialized-const-reference
  CXXFLAGS_WARNINGS += -Wunreachable-code
  CXXFLAGS_WARNINGS += -Wunused
  CXXFLAGS_WARNINGS += -Wzero-as-null-pointer-constant
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 
  CXXFLAGS_WARNINGS += 

  #CXXFLAGS_WARNINGS += -Wno-unused-variable


	#Some of the .... suchen auf der CLANG seite 

	# TODO: Die unused lambda captures will ersetzen durch generelle camptures 

endif
 
endif
 
CXXFLAGS_WARNINGS += -Wno-conversion -Wno-sign-compare -Wno-unused-variable -Wno-unused-parameter -Wno-vla -Wno-type-limits 




 
 
 
 
 
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

SANITIZER_FLAG := -pg -fno-omit-frame-pointer -ftrapv

ifeq ($(FLAG_DO_USE_SANITIZER),yes)

	SANITIZERS :=

	ifeq ($(FLAG_CXX),GCC)

		SANITIZERS :=$(SANITIZERS)pointer-compare,pointer-subtract,
		SANITIZERS :=$(SANITIZERS)undefined,

	else ifeq ($(FLAG_CXX),CLANG)

		SANITIZERS :=$(SANITIZERS)pointer-compare,pointer-subtract,
		SANITIZERS :=$(SANITIZERS)undefined,

		SANITIZERS :=$(SANITIZERS)float-divide-by-zero,
		SANITIZERS :=$(SANITIZERS)unsigned-integer-overflow,
		SANITIZERS :=$(SANITIZERS)implicit-conversion,
		SANITIZERS :=$(SANITIZERS)nullability-arg,
		SANITIZERS :=$(SANITIZERS)nullability-assign,
		SANITIZERS :=$(SANITIZERS)nullability-return,

		#SANITIZERS :=$(SANITIZERS)memory,

	endif

	SANITIZERS :=$(SANITIZERS)address,leak,
	# SANITIZERS :=$(SANITIZERS)thread 

	# thread cannot be combined with address and leak 

	
	# Comment out the following line to disable ALL built-in sanitizers 
	SANITIZER_FLAG += -fsanitize=$(SANITIZERS) 

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

CPPFLAGS := $(FLAG_DO_NOT_CHECK_MESHES) $(FLAG_DONTASSERT) $(FLAG_DO_COMPILEDEBUGMODE)



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




