
# 
# This file contains definitions of variables 
# to be used in the compilation process 
# for objects and shared libraries 
# 
# This file produces the variables 
# 	- CXX
# 	- CXXFLAGS
# 	- CPPFLAGS
# 	- LDLIBS
# 
# Note that these include several preprocessor variables
# Furthermore, it sets:
# 
#   - LINKINGTYPE
# 
#  






# Do you want to use GCC or Clang?
# Uncomment the appropriate definition below
# FLAG_CXX := CLANG
FLAG_CXX := GCC
# FLAG_CXX := ICC


# Do you want to DISABLE the general assert macro?
# Uncomment the following line to disable the general assert macro
# FLAG_DISABLE_ASSERTIONS=yes

# Do you want the standard library assert macro instead of the custom one?
# Uncomment the following line to use the standard library assert macro 
# FLAG_USE_ORIGINAL_ASSERT_MACRO=yes

# Do you want to DISABLE checking of meshes?
# Uncomment the following line to disable extensive check routines for meshes
FLAG_DISABLE_CHECK_MESHES=yes

# Do you want to DISABLE the custom logging framework
# in favor of standard library routines?
# Uncomment the following line for that
# FLAG_USE_PRIMITIVE_LOGGING=yes

# Do you want to ENABLE the standard library debugging flags 
# Uncomment the following line to enable the standard library debugging flags 
# FLAG_DISABLE_STDLIBDEBUG=yes

# Do you want to DISABLE excpetion handling?
# Uncomment the following line to disable exception handling
# FLAG_NO_EXCEPTIONS=yes

# Do you want to ENABLE extended precision?
# Uncomment the following line to switch from double precision to extended precision
# FLAG_DO_USE_EXTENDED_PRECISION=yes 
# := -DEXTENDED_PRECISION

# Do you want to ENABLE excessive warning options?
# Uncomment the following line to enable excessive warning options
# FLAG_EXCESSIVE_WARNINGS=yes

# Do you want to ENABLE the Clang sanitizer?
# Uncomment the following line to enable compilation with the Clang sanitizer
# FLAG_DO_USE_SANITIZER=yes

# Do you want to enable static analysis during the compilation process
# Uncomment the following line to enable static analysis
# FLAG_DO_STATICANALYSIS=yes

# Do you want to compile with all optimization flags enabled?
# Uncomment the following line to have this done so
# FLAG_DO_OPTIMIZE=yes

# Do you want to ENABLE the use of openMP?
# Uncomment the following line to enable compilation with openMP
# FLAG_ENABLE_OPENMP=yes

# Do you want to ENABLE the use of tcmalloc?
# Uncomment the following line to enable tcmalloc
# FLAG_USE_TCMALLOC=yes

# Do you want to DISABLE embedding of Debug information?
# Uncomment the following line to have no debug information included
# FLAG_NO_DEBUGINFO=yes

# Do you want to ENABLE the backtracer for debugging?
# Uncomment the following line to enable backtracing
# FLAG_USE_BACKTRACER=yes

# Do you want to strip unused symbols from the executables?
# Uncomment the following line to accomplish this
# FLAG_DO_STRIP=yes

# Do you want to ENABLE profile generation? 
# Uncomment the following line to enable profile generation at every run.
# FLAG_DO_PROFILE=yes





# If we are in RELEASE_MODE then set the following flags 

ifdef $(RELEASE_MODE)
FLAG_DISABLE_CHECK_MESHES=yes
FLAG_DISABLE_STDLIBDEBUG=yes
FLAG_DISABLE_ASSERTIONS=yes
FLAG_DO_OPTIMIZE=yes
endif




# print all the compilation flags set manually or automatically
.PHONY: parameters
parameters:
	$(info FLAG_CXX                       = $(FLAG_CXX) )
	$(info FLAG_DISABLE_ASSERTIONS        = $(FLAG_DISABLE_ASSERTIONS) ) 
	$(info FLAG_USE_ORIGINAL_ASSERT_MACRO = $(FLAG_USE_ORIGINAL_ASSERT_MACRO) )
	$(info FLAG_DISABLE_CHECK_MESHES      = $(FLAG_DISABLE_CHECK_MESHES) ) 
	$(info FLAG_USE_PRIMITIVE_LOGGING     = $(FLAG_USE_PRIMITIVE_LOGGING) )
	$(info FLAG_DISABLE_STDLIBDEBUG       = $(FLAG_DISABLE_STDLIBDEBUG) ) 
	$(info FLAG_NO_EXCEPTIONS             = $(FLAG_NO_EXCEPTIONS) ) 
	$(info FLAG_DO_USE_EXTENDED_PRECISION = $(FLAG_DO_USE_EXTENDED_PRECISION) ) 
	$(info FLAG_EXCESSIVE_WARNINGS        = $(FLAG_EXCESSIVE_WARNINGS) ) 
	$(info FLAG_DO_USE_SANITIZER          = $(FLAG_DO_USE_SANITIZER) ) 
	$(info FLAG_DO_STATICANALYSIS         = $(FLAG_DO_STATICANALYSIS) ) 
	$(info FLAG_DO_OPTIMIZE               = $(FLAG_DO_OPTIMIZE) ) 
	$(info FLAG_ENABLE_OPENMP             = $(FLAG_ENABLE_OPENMP) ) 
	$(info FLAG_USE_TCMALLOC              = $(FLAG_USE_TCMALLOC) ) 
	$(info FLAG_NO_DEBUGINFO              = $(FLAG_NO_DEBUGINFO) ) 
	$(info FLAG_USE_BACKTRACER            = $(FLAG_USE_BACKTRACER) ) 
	$(info FLAG_DO_STRIP                  = $(FLAG_DO_STRIP) ) 
	$(info FLAG_DO_PROFILE                = $(FLAG_DO_PROFILE) ) 
	@true










###############################################
#                                             #
#         Set the compiler command            #
#                                             #
###############################################

ifeq ($(FLAG_CXX),GCC)

  CXX := g++ -std=c++2a
  
else ifeq ($(FLAG_CXX),CLANG)

  CXX := clang++ -std=c++2a

else ifeq ($(FLAG_CXX),ICC)

  CXX := icc -std=c++2a

else

  $(error No compiler recognized)

endif








###############################################
#                                             #
#           Language std settings             #
#                                             #
###############################################

CXXFLAGS_LANG := -std=c++2a -pedantic 












###############################################
#                                             #
#               Optimization                  #
#                                             #
###############################################

# If optimization is enabled, then set a number of flags 
# In the absence of optimization, we set O1

CXXFLAGS_OPTIMIZE:=

ifeq ($(FLAG_DO_OPTIMIZE),yes)

	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_OPTIMIZE += -march=core-avx2
		CXXFLAGS_OPTIMIZE += -intel-optimized-headers 
		CXXFLAGS_OPTIMIZE += -Ofast -xHOST -O3 -ipo -no-prec-div -fp-model fast=2
		CXXFLAGS_OPTIMIZE += 
		CXXFLAGS_OPTIMIZE += 
		CXXFLAGS_OPTIMIZE += 
	else 
		CXXFLAGS_OPTIMIZE += -Ofast 
		CXXFLAGS_OPTIMIZE += -march=native 
		ifeq ($(FLAG_CXX),CLANG)
			CXXFLAGS_OPTIMIZE := -inline-threshold=1200
		endif
		CXXFLAGS_OPTIMIZE += -flto
	endif
else
	CXXFLAGS_OPTIMIZE += -O1
endif


# Do we apply OpenMP?

ifeq ($(FLAG_ENABLE_OPENMP),yes)
	CXXFLAGS_OPTIMIZE += -fopenmp
endif


# Do we strip debug information?

ifeq ($(FLAG_CXX),GCC) 
ifeq ($(FLAG_DO_STRIP),yes)
	CXXFLAGS_OPTIMIZE += -ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--strip-all 
endif
endif




###############################################
#                                             #
#        Misc Code generation options         #
#                                             #
###############################################

CXXFLAGS_CODEGEN := 

ifneq ($(FLAG_NO_EXCEPTIONS),yes)
	CXXFLAGS_CODEGEN += -fno-exceptions
endif

CXXFLAGS_CODEGEN += -fvisibility=default

ifneq ($(OS),Windows_NT)
	CXXFLAGS_CODEGEN += -fpic 
endif



















###############################################
#                                             #
#        Compiler warning settings            #
#                                             #
###############################################

CXXFLAGS_WARNINGS := 
CXXFLAGS_WARNINGS += -Wall -Wextra -Wpedantic 

ifeq ($(FLAG_EXCESSIVE_WARNINGS),yes)

	CXXFLAGS_WARNINGS += -Wodr -Wmissing-field-initializers -Wctor-dtor-privacy -Wsign-promo -Woverloaded-virtual -Wno-missing-braces

	CXXFLAGS_WARNINGS += -Wundef 
	CXXFLAGS_WARNINGS += -Wcast-align -Wcast-qual -Wmissing-declarations -Wredundant-decls -Wno-redundant-decls
	CXXFLAGS_WARNINGS += -Wformat=2 -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wmisleading-indentation
	# CXXFLAGS_WARNINGS += -Wold-style-cast

	ifeq ($(FLAG_CXX),GCC)

		CXXFLAGS_WARNINGS += -Wmultiple-inheritance -Wvirtual-inheritance
		CXXFLAGS_WARNINGS += -Wpointer-arith
		CXXFLAGS_WARNINGS += -Wreturn-local-addr
		CXXFLAGS_WARNINGS += -Wfloat-equal -Wdouble-promotion -Wfloat-conversion
		CXXFLAGS_WARNINGS += -Wno-sign-compare
		CXXFLAGS_WARNINGS += -Wconversion 
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
		CXXFLAGS_WARNINGS += -Wunsafe-loop-optimizations 
		#   CXXFLAGS_WARNINGS += -Winline 
		CXXFLAGS_WARNINGS += -Wvector-operation-performance
		CXXFLAGS_WARNINGS += -Wdisabled-optimization
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
		#   CXXFLAGS_WARNINGS += -Wdouble-promotion #disabled
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
		CXXFLAGS_WARNINGS += -Wreserved-id-macro
		CXXFLAGS_WARNINGS += -Wreserved-user-defined-literal
		CXXFLAGS_WARNINGS += -Wreturn-std-move
		CXXFLAGS_WARNINGS += -Wself-assign -Wself-move
		CXXFLAGS_WARNINGS += -Wsemicolon-before-method-body
		CXXFLAGS_WARNINGS += -Wsometimes-uninitialized
		CXXFLAGS_WARNINGS += -Wstrict-prototypes
		CXXFLAGS_WARNINGS += -Wstring-conversion
		#   CXXFLAGS_WARNINGS += -Wsuggest-destructor-override #unknown
		CXXFLAGS_WARNINGS += -Wtautological-compare
		CXXFLAGS_WARNINGS += -Wtautological-type-limit-compare
		#   CXXFLAGS_WARNINGS += -Wuninitialized-const-reference #unknown
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


# In any case, remove the following warnings 

CXXFLAGS_WARNINGS += -Wno-conversion 
CXXFLAGS_WARNINGS += -Wno-sign-compare 
CXXFLAGS_WARNINGS += -Wno-unused-variable
CXXFLAGS_WARNINGS += -Wno-unused-parameter
CXXFLAGS_WARNINGS += -Wno-vla
CXXFLAGS_WARNINGS += -Wno-unknown-pragmas
CXXFLAGS_WARNINGS += -Wno-type-limits 















###############################################
#                                             #
#          Debugging information              #
#                                             #
###############################################

CXXFLAGS_DEBUG :=
ifeq ($(FLAG_NO_DEBUGINFO),yes)
else
CXXFLAGS_DEBUG += -g
endif


###############################################
#                                             #
#          Profiling instrumentation          #
#                                             #
###############################################

ifeq ($(FLAG_DO_PROFILE),yes)

	CXXFLAGS_PROF:= -pg -fno-omit-frame-pointer 

endif


###############################################
#                                             #
#          Sanitizer instrumentation          #
#                                             #
###############################################

# There are several incompatibilities between sanitizers 
# thread cannot be combined with address and leak
# Warning: memory causes considerable slowdown

ifeq ($(FLAG_DO_USE_SANITIZER),yes)

	SANITIZERS := -ftrapv 

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
	CXXFLAGS_SANI += -fsanitize=$(SANITIZERS) 

endif


###############################################
#                                             #
#    Use TCMalloc instead of std allocators   #
#                                             #
###############################################

ifeq ($(FLAG_USE_TCMALLOC),yes)
	CXXFLAGS_MALLOC=-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free
else
	CXXFLAGS_MALLOC=
endif








 
 
###############################################
#                                             #
#          Static analysis                    #
#                                             #
###############################################

ifeq ($(FLAG_DO_STATICANALYSIS),yes)

	ifeq ($(FLAG_CXX),GCC)

		CXXFLAGS_STATICANALYSER := -fanalyzer -Wanalyzer-too-complex

	else ifeq ($(FLAG_CXX),CLANG)

		CXXFLAGS_STATICANALYSER := 

	endif

endif


###############################################
#                                             #
#        Format of diagnostic settings        #
#                                             #
###############################################

ifeq ($(FLAG_CXX),CLANG)

  CXXFLAGS_DIAGNOSISFORMAT := -fdiagnostics-show-template-tree

endif

















###############################################
#                                             #
#   CXXFLAGS - gather C++ compiler flags      #
#                                             #
###############################################

CXXFLAGS := 
CXXFLAGS += ${CXXFLAGS_LANG}
CXXFLAGS += ${CXXFLAGS_DIAGNOSISFORMAT}
CXXFLAGS += ${CXXFLAGS_WARNINGS}
CXXFLAGS += ${CXXFLAGS_STATICANALYSER}
CXXFLAGS += ${CXXFLAGS_DEBUG}
CXXFLAGS += $(CXXFLAGS_PROF)
CXXFLAGS += $(CXXFLAGS_SANI)
CXXFLAGS += ${CXXFLAGS_MALLOC}
CXXFLAGS += ${CXXFLAGS_OPTIMIZE}
CXXFLAGS += ${CXXFLAGS_CODEGEN}

CXXFLAGS := $(strip $(CXXFLAGS))


CXXFLAGS_EXECUTABLE:=
CXXFLAGS_EXECUTABLE+=$(CXXFLAGS)

ifeq ($(FLAG_DO_OPTIMIZE),yes)
	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_EXECUTABLE+=-fwhole-program
	else 
		ifeq ($(FLAG_CXX),GCC)
			CXXFLAGS_EXECUTABLE+=-fwhole-program
		endif
	endif
else
endif



###############################################
#                                             #
#   CPPFLAGS - gather preprocessor flags      #
#                                             #
###############################################

CPPFLAGS := 

ifeq ($(FLAG_DISABLE_CHECK_MESHES),yes)
CPPFLAGS += -DDO_NOT_CHECK_MESHES
endif

ifeq ($(FLAG_USE_BACKTRACER),yes)
CPPFLAGS += -DUSE_BACKTRACER
endif

ifeq ($(FLAG_DISABLE_ASSERTIONS),yes)
CPPFLAGS += -DNDEBUG
endif

ifeq ($(FLAG_DISABLE_STDLIBDEBUG),yes)
else 
CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_ASSERTIONS -D_GLIBCXX_SANITIZE_VECTOR
endif

ifeq ($(FLAG_DO_USE_EXTENDED_PRECISION),yes)
CPPFLAGS += -DEXTENDED_PRECISION
endif

CPPFLAGS := $(strip $(CPPFLAGS))


###############################################
#                                             #
#         LDLIBS - gather linker flags        #
#                                             #
###############################################

ifeq ($(FLAG_USE_TCMALLOC),yes)
	LDLIBS :=-l:libtcmalloc.so.4
else
	LDLIBS :=
endif

LDLIBS := $(strip $(LDLIBS))






###############################################
#                                             #
#         Choose the type of linking          #
#         - static                            #
#         - dynamic                           #
#         - unspecified                       #
#                                             #
###############################################


LINKINGTYPE:=unspecified

ifeq ($(OS),Windows_NT)
LINKINGTYPE:=static
else ifeq ($(FLAG_DO_OPTIMIZE=yes),yes)
LINKINGTYPE:=static
else
LINKINGTYPE:=dynamic
endif




