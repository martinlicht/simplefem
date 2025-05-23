
# 
# This file contains definitions for the compilation
# of objects and shared libraries.
# 
# This file produces the variables:
# 	- CXX
# 	- CXXFLAGS
# 	- CPPFLAGS
# 	- LDFLAGS
# 	- LDLIBS
# 
# These are assembled from several preprocessor variables.
# Furthermore, it sets:
# 
#   - LINKINGTYPE
# 
# 

################################################################################
# EXPECTED VARIABLES:
# - projectdir : path/to/project/directory
ifndef projectdir
$(error Expect 'projectdir')
endif





# Uncomment your choice of compiler below
FLAG_CXX := CLANG
# FLAG_CXX := GCC
# FLAG_CXX := ICC


# Do you want to DISABLE the general assert macro?
# Uncomment the following line to disable the general assert macro:
# FLAG_DISABLE_ASSERTIONS=yes

# Do you want the standard library assert macro instead of the custom one?
# Uncomment the following line to use the standard library assert macro:
# FLAG_USE_ORIGINAL_ASSERT_MACRO=yes

# Do you want assert messages to be discarded?
# Uncomment the following line to simplify the debugging macros:
# FLAG_DISCARD_ASSERT_MESSAGES=yes

# Do you want to DISABLE checking of meshes?
# Uncomment the following line to disable extensive check routines for meshes:
FLAG_DISABLE_CHECK_MESHES=yes

# Do you want to ENABLE the standard library debugging flags?
# Uncomment the following line to enable the standard library debugging flags:
# FLAG_DISABLE_STDLIBDEBUG=yes

# Do you want to DISABLE the custom logging framework
# in favor of standard library routines?
# Uncomment the following line for that:
# FLAG_USE_PRIMITIVE_LOGGING=yes

# Do you want to DISABLE excpetion handling?
# Uncomment the following line to disable exception handling:
FLAG_NO_EXCEPTIONS=yes

# Do you want to compile with all optimization flags enabled?
# Uncomment the following line to have this done so:
FLAG_DO_OPTIMIZE=yes

# Do you want to ENABLE the use of openMP?
# Uncomment the following line to enable compilation with openMP:
FLAG_ENABLE_OPENMP=yes

# Do you want to ENABLE excessive warning options?
# Uncomment the following line to enable excessive warning options:
FLAG_EXCESSIVE_WARNINGS=yes

# Do you want to ENABLE either extended precision or single precision?
# Uncomment one the following lines to switch from double precision
# to either extended precision or single precision:
# FLAG_DO_USE_EXTENDED_PRECISION=yes
# FLAG_DO_USE_SINGLE_PRECISION=yes

# Logging output in color?
# FLAG_COLORED_OUTPUT=yes

# Do you want to ENABLE the Clang sanitizer?
# Uncomment the following line to enable compilation with the Clang sanitizer:
# FLAG_DO_USE_SANITIZER=yes

# Do you want to enable static analysis during the compilation process?
# Uncomment the following line to enable static analysis:
# FLAG_DO_STATICANALYSIS=yes

# Do you want to ENABLE the use of TCMalloc?
# Uncomment the following line to enable TCMalloc:
# FLAG_USE_TCMALLOC=yes

# Do you want to DISABLE embedding of Debug information?
# Uncomment the following line to have no debug information included:
# FLAG_NO_DEBUGINFO=yes

# Do you want to ENABLE the backtracer for debugging (available since C++23)?
# Uncomment the following line to enable backtracing:
# FLAG_USE_BACKTRACER=yes

# Do you want to strip unused symbols from the executables?
# Uncomment the following line to accomplish this:
# FLAG_DO_STRIP=yes

# Do you want to ENABLE profile generation?
# Uncomment the following line to enable profile generation at every run:
# FLAG_DO_PROFILE=yes

# Choose the linker by uncommenting, or leave commented for default linker:
# LDFLAGS += -fuse-ld=bfd
# LDFLAGS += -fuse-ld=lld
# LDFLAGS += -fuse-ld=gold
# LDFLAGS += -fuse-ld=mold





###############################################
#                                             #
#       Post-process the options              #
#                                             #
###############################################

# Use this file to overwrite the default settings above on a local machine.
# Here, this will only affect whether any of the flags above is set/unset,
# since any strings will be overwritten later in this file.
-include $(projectdir)/OVERWRITE.COMPILE.mk

# If we are in RELEASE_MODE (variable set on commandline or elsewhere),
# then set the following flags.

ifdef OPTIMIZED_DEBUG
#FLAG_DISABLE_CHECK_MESHES=yes
#FLAG_DISABLE_STDLIBDEBUG=yes
#FLAG_DISABLE_ASSERTIONS=yes
FLAG_DO_OPTIMIZE=yes
FLAG_ENABLE_OPENMP=yes
FLAG_NO_EXCEPTIONS=yes
endif

ifdef RELEASE_MODE
FLAG_DISABLE_CHECK_MESHES=yes
FLAG_DISABLE_STDLIBDEBUG=yes
FLAG_DISABLE_ASSERTIONS=yes
FLAG_DO_OPTIMIZE=yes
FLAG_ENABLE_OPENMP=yes
FLAG_NO_DEBUGINFO=yes
FLAG_NO_EXCEPTIONS=yes
FLAG_DO_STRIP=yes
FLAG_USE_PRIMITIVE_LOGGING=yes
endif

# You cannot request single and extended precision at the same time.
ifdef FLAG_DO_USE_EXTENDED_PRECISION
ifdef FLAG_DO_USE_SINGLE_PRECISION
$(error Extended and single precision requested at the same time)
endif
endif

# Single-precision is only supported on GCC.
ifneq ($(FLAG_CXX),GCC)
ifdef FLAG_DO_USE_SINGLE_PRECISION
$(error Single precision floating-point setup only with GCC)
endif
endif





###############################################
#                                             #
#    Set the compiler and the std library     #
#       (see also language std below)         #
#                                             #
###############################################

ifeq ($(FLAG_CXX),GCC)

  CXX := g++
  #-ftime-report

else ifeq ($(FLAG_CXX),CLANG)

  CXX := clang++ -stdlib=libstdc++ 
  # -ftime-trace

else ifeq ($(FLAG_CXX),ICC)

  CXX := icc 

else

  $(error No compiler recognized)

endif

###############################################
#                                             #
#           Language std settings             #
#                                             #
###############################################

CXXFLAGS_LANG := -std=c++14 -pedantic -fno-rtti 

ifeq ($(FLAG_NO_EXCEPTIONS),yes)
	CXXFLAGS_LANG += -fno-exceptions 
endif

# The following optimizations may or may not be beneficial.
# -fno-omit-frame-pointer # The frame pointer is necessary for profiling.
# -fno-unwind-tables	  # Unwind tables are necessary for exceptions.

# If single precision is requested, then all constants become `float` by default.
ifeq ($(FLAG_DO_USE_SINGLE_PRECISION),yes)
	CXXFLAGS_LANG += -fsingle-precision-constant
endif

# If at some point we want to play with GCC's visibility feature.
# CXXFLAGS_LANG += -fvisibility=default

# Inlined class member functions cannot be referenced outside of their translation unit.
# This is supposed to be reduce some of the clutter.
CXXFLAGS_LANG += -fvisibility-inlines-hidden

# enums can only take values in the enumeration type, everything else is undefined behavior
CXXFLAGS_LANG += -fstrict-enums

# Enforce a strict evaluation order. This settles some implementation-defined behavior.
# CXXFLAGS_LANG += -fstrong-eval-order

# Do we enable OpenMP?
ifeq ($(FLAG_ENABLE_OPENMP),yes)
	CXXFLAGS_LANG += -fopenmp
endif










###############################################
#                                             #
#               Optimization                  #
#                                             #
###############################################

# If optimization is enabled, then we set the following flags.
# In the absence of optimization, we set: O0
# 
# Irrespective of the compiler, the following are save to use:
# -fshort-enums: 		 use short enums
# -fno-finite-math-only: do NOT assume that math is finite

CXXFLAGS_OPTIMIZE:=

ifeq ($(FLAG_DO_OPTIMIZE),yes)

	CXXFLAGS_OPTIMIZE += -march=native -mtune=native 

	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_OPTIMIZE += -march=core-avx2
		CXXFLAGS_OPTIMIZE += -intel-optimized-headers 
		CXXFLAGS_OPTIMIZE += -xHOST -O3 -Ofast -ipo -no-prec-div -fp-model fast=2
	else 
# wierd warnings appear at LTO and O1+ ....
		CXXFLAGS_OPTIMIZE += -flto
		CXXFLAGS_OPTIMIZE += -O3
		CXXFLAGS_OPTIMIZE += -fshort-enums
		ifneq ($(FLAG_DO_PROFILE),yes)
			CXXFLAGS_OPTIMIZE += -fomit-frame-pointer
		endif
		ifeq ($(FLAG_CXX),GCC)
#			CXXFLAGS_OPTIMIZE += -finline-limit=1200
			CXXFLAGS_OPTIMIZE += -fno-fat-lto-objects
			CXXFLAGS_OPTIMIZE += -fdevirtualize-speculatively
# 			CXXFLAGS_OPTIMIZE += -fno-signed-zeros -fno-signaling-nans -fno-trapping-math -fassociative-math -fno-math-errno  -funsafe-math-optimizations 
# 			CXXFLAGS_OPTIMIZE += -fmerge-all-constants
# 			CXXFLAGS_OPTIMIZE += -fipa-pta 
# Loop unrolling is very speculative			
			CXXFLAGS_OPTIMIZE += -funroll-loops -fvariable-expansion-in-unroller -floop-nest-optimize
#			CXXFLAGS_OPTIMIZE += -malign-double 
		endif
		ifeq ($(FLAG_CXX),CLANG)
#			 CXXFLAGS_OPTIMIZE += -finline-limit=1200
		endif
	endif
	 
else
	CXXFLAGS_OPTIMIZE += -O0
endif

CXXFLAGS_OPTIMIZE += -fshort-enums
CXXFLAGS_OPTIMIZE += -fno-finite-math-only







###############################################
#                                             #
#        Misc Code generation options         #
#                                             #
###############################################

CXXFLAGS_CODEGEN := 

# If we are NOT on Windows, then:
# -fpic:    use position-independent code,
# -fno-plt: avoid the procedure linkage table.
ifneq ($(OS),Windows_NT)
	CXXFLAGS_CODEGEN += -fpic -fno-plt 
endif






















###############################################
#                                             #
#        Compiler warning settings            #
#                                             #
###############################################

# We generally enable the following warnings.

CXXFLAGS_WARNINGS := 
CXXFLAGS_WARNINGS += -Wall -Wextra -Wpedantic 

CXXFLAGS_WARNINGS += -Wformat=2
CXXFLAGS_WARNINGS += -Wformat-nonliteral
CXXFLAGS_WARNINGS += -Wformat-security
CXXFLAGS_WARNINGS += -Wformat-y2k

	
ifeq ($(FLAG_EXCESSIVE_WARNINGS),yes)

# If excessive warnings are requested, we distinguish between compilers.
# We first enable as much as we can, and then disable anything we do not want.

	ifeq ($(FLAG_CXX),GCC)

		# This flag is buggy and should not be used
		# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=85783
		# CXXFLAGS_WARNINGS += -Wno-alloc-size-larger-than
		
		# # DISABLED: CXXFLAGS_WARNINGS += -Wabi
		CXXFLAGS_WARNINGS += -Waggressive-loop-optimizations
		CXXFLAGS_WARNINGS += -Walloca
		# # DISABLED: CXXFLAGS_WARNINGS += -Walloc-size
		CXXFLAGS_WARNINGS += -Walloc-zero
		CXXFLAGS_WARNINGS += -Warith-conversion
		CXXFLAGS_WARNINGS += -Warray-bounds
		CXXFLAGS_WARNINGS += -Warray-bounds=2
		CXXFLAGS_WARNINGS += -Wattribute-alias=2
		# # DISABLED: CXXFLAGS_WARNINGS += -Wattribute-warning
		CXXFLAGS_WARNINGS += -Wno-attribute-warning
		
		# # unknown: CXXFLAGS_WARNINGS += -Wbidi-chars=any
		
		# # unknown: CXXFLAGS_WARNINGS += -Wcalloc-transposed-args
		CXXFLAGS_WARNINGS += -Wcast-qual
		CXXFLAGS_WARNINGS += -Wcast-align
		CXXFLAGS_WARNINGS += -Wcast-align=strict
		CXXFLAGS_WARNINGS += -Wcomma-subscript
		CXXFLAGS_WARNINGS += -Wconversion
		
		CXXFLAGS_WARNINGS += -Wdangling-else
		CXXFLAGS_WARNINGS += -Wdangling-reference
		CXXFLAGS_WARNINGS += -Wdate-time
		# CXXFLAGS_WARNINGS += -Wdeprecated-variadic-comma-omission
		CXXFLAGS_WARNINGS += -Wdisabled-optimization
		CXXFLAGS_WARNINGS += -Wduplicated-branches
		CXXFLAGS_WARNINGS += -Wduplicated-cond
		
		CXXFLAGS_WARNINGS += -Wexpansion-to-defined
		
		CXXFLAGS_WARNINGS += -Wformat=2
		CXXFLAGS_WARNINGS += -Wformat-overflow=2
		CXXFLAGS_WARNINGS += -Wformat-signedness
		CXXFLAGS_WARNINGS += -Wformat-truncation=1
		# # disabled: CXXFLAGS_WARNINGS += -Wformat-truncation=2
		CXXFLAGS_WARNINGS += -Wfloat-conversion
		CXXFLAGS_WARNINGS += -Wfloat-equal
		CXXFLAGS_WARNINGS += -Wfree-nonheap-object
		
		CXXFLAGS_WARNINGS += -Wimplicit-fallthrough=4
		# # DISABLED: CXXFLAGS_WARNINGS += -Winfinite-recursion
		CXXFLAGS_WARNINGS += -Winit-self
		# # DISABLED: CXXFLAGS_WARNINGS += -Winline
		# # DISABLED: CXXFLAGS_WARNINGS += -Wno-inline
		CXXFLAGS_WARNINGS += -Wint-in-bool-context
		CXXFLAGS_WARNINGS += -Winvalid-pch
		# # unknown: CXXFLAGS_WARNINGS += -Winvalid-utf8
		
		CXXFLAGS_WARNINGS += -Wlogical-op
		# # DISABLED: CXXFLAGS_WARNINGS += -Wlong-long
		
		# # DISABLED: CXXFLAGS_WARNINGS += -Wmisleading-indentation
		CXXFLAGS_WARNINGS += -Wmissing-declarations
		CXXFLAGS_WARNINGS += -Wmissing-include-dirs
		CXXFLAGS_WARNINGS += -Wmultiple-inheritance
		
		# # DISABLED: CXXFLAGS_WARNINGS += -Wnrvo
		CXXFLAGS_WARNINGS += -Wnull-dereference
		
		CXXFLAGS_WARNINGS += -Woverloaded-virtual=2
		
		# # DISABLED: CXXFLAGS_WARNINGS += -Wpadded
		CXXFLAGS_WARNINGS += -Wparentheses
		CXXFLAGS_WARNINGS += -Wpointer-arith
		CXXFLAGS_WARNINGS += -Wpointer-arith
		
		CXXFLAGS_WARNINGS += -Wredundant-decls
		CXXFLAGS_WARNINGS += -Wreturn-local-addr
		# # DISABLED: CXXFLAGS_WARNINGS += -Wreturn-mismatch
		
		# # DISABLED: CXXFLAGS_WARNINGS += -Wshadow
		# # DISABLED: CXXFLAGS_WARNINGS += -Wshadow=global
		CXXFLAGS_WARNINGS += -Wshift-negative-value
		CXXFLAGS_WARNINGS += -Wshift-overflow=2
		CXXFLAGS_WARNINGS += -Wsign-compare
		CXXFLAGS_WARNINGS += -Wsign-conversion
		CXXFLAGS_WARNINGS += -Wsign-promo
		CXXFLAGS_WARNINGS += -Wstrict-aliasing=3
		# # DISABLED: CXXFLAGS_WARNINGS += -Wstrict-flex-arrays
		# # DISABLED: CXXFLAGS_WARNINGS += -Wstrict-overflow=5
		CXXFLAGS_WARNINGS += -Wstringop-overflow=4
		# # DISABLED: CXXFLAGS_WARNINGS += -Wsuggest-attribute
		# # DISABLED: CXXFLAGS_WARNINGS += -Wsuggest-final-methods
		CXXFLAGS_WARNINGS += -Wno-suggest-final-methods
		# # DISABLED: CXXFLAGS_WARNINGS += -Wsuggest-final-types
		CXXFLAGS_WARNINGS += -Wno-suggest-final-types
		CXXFLAGS_WARNINGS += -Wsuggest-override
		CXXFLAGS_WARNINGS += -Wswitch-default
		CXXFLAGS_WARNINGS += -Wswitch-enum
		CXXFLAGS_WARNINGS += -Wsync-nand
		
		CXXFLAGS_WARNINGS += -Wtautological-compare
		# # DISABELD: CXXFLAGS_WARNINGS += -Wtraditional
		CXXFLAGS_WARNINGS += -Wtrampolines
		# # DISABELD: CXXFLAGS_WARNINGS += -Wtrivial-auto-var-init
		
		CXXFLAGS_WARNINGS += -Wuninitialized
		CXXFLAGS_WARNINGS += -Wunused
		CXXFLAGS_WARNINGS += -Wunused-const-variable=1
		# # DISABLED: CXXFLAGS_WARNINGS += -Wunused-macros
		# # DISABLED: CXXFLAGS_WARNINGS += -Wuseless-cast
		CXXFLAGS_WARNINGS += -Wunsafe-loop-optimizations
		CXXFLAGS_WARNINGS += -Wundef
		
		CXXFLAGS_WARNINGS += -Wvector-operation-performance
		CXXFLAGS_WARNINGS += -Wvirtual-inheritance
		
		CXXFLAGS_WARNINGS += -Wwrite-strings
		
		CXXFLAGS_WARNINGS += -Wzero-length-bounds
		CXXFLAGS_WARNINGS += -Wzero-as-null-pointer-constant
		
		# # check which options there are ...
		# # TODO What is stack smashing???? HSA????
		# # TODO Read the format warnings

		# CXXFLAGS_WARNINGS += -Wshadow			# we accept shadowed variables
		CXXFLAGS_WARNINGS += -Wshadow-local			# ... except local variables shadowing other locals


		CXXFLAGS_WARNINGS += 
		
	else ifeq ($(FLAG_CXX),CLANG)

		CXXFLAGS_WARNINGS += -Weverything
		
		# Clang only
		CXXFLAGS_WARNINGS += -Wno-binary-literal
		CXXFLAGS_WARNINGS += -Wno-cast-function-type-strict
		CXXFLAGS_WARNINGS += -Wno-c++98-compat
		CXXFLAGS_WARNINGS += -Wno-c++98-compat-bind-to-temporary-copy
		CXXFLAGS_WARNINGS += -Wno-c++98-compat-pedantic
		CXXFLAGS_WARNINGS += -Wno-c++-compat
		CXXFLAGS_WARNINGS += -Wno-covered-switch-default
		CXXFLAGS_WARNINGS += -Wno-documentation
		CXXFLAGS_WARNINGS += -Wno-double-promotion #disabled
		CXXFLAGS_WARNINGS += -Wno-dtor-name
		CXXFLAGS_WARNINGS += -Wno-extra-semi-stmt
		CXXFLAGS_WARNINGS += -Wno-exit-time-destructors
		CXXFLAGS_WARNINGS += -Wno-float-equal
		CXXFLAGS_WARNINGS += -Wno-frame-larger-than
		CXXFLAGS_WARNINGS += -Wno-global-constructors
		CXXFLAGS_WARNINGS += -Wno-inconsistent-missing-destructor-override
		CXXFLAGS_WARNINGS += -Wno-invalid-utf8
		CXXFLAGS_WARNINGS += -Wno-local-type-template-args
		CXXFLAGS_WARNINGS += -Wno-misleading-indentation
		CXXFLAGS_WARNINGS += -Wno-missing-prototypes
		CXXFLAGS_WARNINGS += -Wno-old-style-cast
		CXXFLAGS_WARNINGS += -Wno-padded
		CXXFLAGS_WARNINGS += -Wno-shadow
		CXXFLAGS_WARNINGS += -Wno-shadow-all
		CXXFLAGS_WARNINGS += -Wno-shadow-field-in-constructor
		CXXFLAGS_WARNINGS += -Wno-source-uses-openacc
		CXXFLAGS_WARNINGS += -Wno-suggest-destructor-override
		CXXFLAGS_WARNINGS += -Wno-tautological-negation-compare
		CXXFLAGS_WARNINGS += -Wno-uninitialized-const-reference
		CXXFLAGS_WARNINGS += -Wno-unsafe-buffer-usage
		CXXFLAGS_WARNINGS += -Wno-unused-template
		CXXFLAGS_WARNINGS += -Wno-unreachable-code-loop-increment
		CXXFLAGS_WARNINGS += -Wno-weak-vtables
		CXXFLAGS_WARNINGS += 
		CXXFLAGS_WARNINGS += 
		CXXFLAGS_WARNINGS += 
		
		# TODO:(martin) shadowed and unused variables should be revisited
		CXXFLAGS_WARNINGS += -Wno-shadow
		CXXFLAGS_WARNINGS += -Wno-unused-variable
		CXXFLAGS_WARNINGS += -Wno-shadow-field
		CXXFLAGS_WARNINGS += -Wno-unreachable-code
		
		# TODO:(martin) should virtual destructors be marked override?
		CXXFLAGS_WARNINGS += -Wno-inconsistent-missing-destructor-override
		CXXFLAGS_WARNINGS += -Wno-suggest-destructor-override
		CXXFLAGS_WARNINGS += -Wno-documentation-unknown-command
		
		# TODO:(martin) learn about weak vtables
		CXXFLAGS_WARNINGS += -Wno-weak-vtables
		
	endif

endif

# Having set the warnings specific to any compiler,
# the following warnings are enabled for all compilers.

CXXFLAGS_WARNINGS += -Weffc++
CXXFLAGS_WARNINGS += -Wctor-dtor-privacy
CXXFLAGS_WARNINGS += -Wmisleading-indentation
CXXFLAGS_WARNINGS += -Wmissing-declarations
CXXFLAGS_WARNINGS += -Wmissing-field-initializers
CXXFLAGS_WARNINGS += -Wnon-virtual-dtor
CXXFLAGS_WARNINGS += -Woverloaded-virtual
CXXFLAGS_WARNINGS += -Wodr
CXXFLAGS_WARNINGS += -Wredundant-decls
CXXFLAGS_WARNINGS += -Wsign-promo
CXXFLAGS_WARNINGS += -Wundef
CXXFLAGS_WARNINGS += -Wunused

CXXFLAGS_WARNINGS += -Wno-old-style-cast

# possible experimental feature: CXXFLAGS_WARNINGS += -Wlifetime
# https://stackoverflow.com/questions/52662135/the-purpose-of-the-wlifetime-flag


# Finally, the following warnings are disabled in any case.

CXXFLAGS_WARNINGS += -Wno-double-promotion      # it is okay to implicitly promote float to double
CXXFLAGS_WARNINGS += -Wno-old-style-cast        # we accept old C-style casts
CXXFLAGS_WARNINGS += -Wno-unused-variable       # we accept unused variables
CXXFLAGS_WARNINGS += -Wno-unused-parameter      # we accept unused parameters
CXXFLAGS_WARNINGS += -Wno-unknown-pragmas       # we accept unknown pragmas
CXXFLAGS_WARNINGS += -Wno-vla                   # we accept variable-length arrays

CXXFLAGS_WARNINGS += -Wno-conversion            # mostly unimportant stuff
CXXFLAGS_WARNINGS += -Wno-float-equal           # numerous comparisons of float to zero
CXXFLAGS_WARNINGS += -Wno-missing-braces        # whether braces are missing anywhere
CXXFLAGS_WARNINGS += -Wno-sign-compare          # numerous comparisons between signed and unsigned indices
CXXFLAGS_WARNINGS += -Wno-sign-conversion       # too many signed indices

# ... and we once again make a compiler distinction.

ifeq ($(FLAG_CXX),GCC)

# CXXFLAGS_WARNINGS += -Wno-stack-usage
CXXFLAGS_WARNINGS += -Wno-c++20-extensions

else ifeq ($(FLAG_CXX),CLANG)

# Defaulted constructors/assignment operators are deleted if the base has them deleted. That's ok.
CXXFLAGS_WARNINGS += -Wno-defaulted-function-deleted

# CXXFLAGS_WARNINGS += -Wno-gnu-zero-variadic-macro-arguments
# CXXFLAGS_WARNINGS += -Wno-vla-extension
# CXXFLAGS_WARNINGS += -Wno-shorten-64-to-32

# We freely use the designated initializer syntax of C++20.
CXXFLAGS_WARNINGS += -Wno-c++20-designator

endif














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

# -pg: add instrumentation for gprof.
# -p:  add instrumentation for prof.
# In addition:
# -fno-omit-frame-pointer: we must keep the frame pointer for profiling.
ifeq ($(FLAG_DO_PROFILE),yes)

	CXXFLAGS_PROF:= -pg -fno-omit-frame-pointer

endif


###############################################
#                                             #
#          Sanitizer instrumentation          #
#                                             #
###############################################

# There are several incompatibilities between sanitizers.
# `thread` cannot be combined with `address` and `leak`.

# NOTE:
# `memory` causes considerable slowdown.

# NOTE:
# The following construction of the strings ensures that there are no spaces added.
# This is necessary for the command line argument.

ifeq ($(FLAG_DO_USE_SANITIZER),yes)

	SANITIZERS :=
	
	ifeq ($(FLAG_CXX),GCC)

		# `address` checks errors such `out-of-bounds` or `use-after-free`
		SANITIZERS :=$(SANITIZERS)address,
		SANITIZERS :=$(SANITIZERS)pointer-compare,
		SANITIZERS :=$(SANITIZERS)pointer-subtract,
		
		# leak detection
		SANITIZERS :=$(SANITIZERS)leak,
		
		# data race detection
		# SANITIZERS :=$(SANITIZERS)thread,

		# check undefined behavior
		SANITIZERS :=$(SANITIZERS)undefined,

	else ifeq ($(FLAG_CXX),CLANG)

		# `address` checks errors such `out-of-bounds` or `use-after-free`
		SANITIZERS :=$(SANITIZERS)address,
		SANITIZERS :=$(SANITIZERS)pointer-compare,
		SANITIZERS :=$(SANITIZERS)pointer-subtract,
		
		# leak detection
		SANITIZERS :=$(SANITIZERS)leak,
		
		# uninitialized memory detector (slowdown ca. 300%)
		SANITIZERS :=$(SANITIZERS)memory,
		
		# data race detection
		# SANITIZERS :=$(SANITIZERS)thread,

		# check undefined behavior
		SANITIZERS :=$(SANITIZERS)undefined,
		
		# The following seem clang-specific and redundant.
		# SANITIZERS :=$(SANITIZERS)float-divide-by-zero,
		# SANITIZERS :=$(SANITIZERS)unsigned-integer-overflow,
		# SANITIZERS :=$(SANITIZERS)implicit-conversion,
		# SANITIZERS :=$(SANITIZERS)nullability-arg,
		# SANITIZERS :=$(SANITIZERS)nullability-assign,
		# SANITIZERS :=$(SANITIZERS)nullability-return,

	endif

	# Detects signed integer over/underflows
	# SANITIZERS :=-ftrapv,

	# Comment out the following line to disable ALL built-in sanitizers
	CXXFLAGS_SANI := -fsanitize=$(SANITIZERS) -ftrapv

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
#        Stripping debug information          #
#                                             #
###############################################

# Do we strip debug information? This may slightly increase compilation time.
# -ffunction-sections: put each functions into a separate section of the binary.
# -fdata-sections:     put each variable into its own section.
# -Wl,--gc-sections:   remove unused sections. Powerful in combination with the previous two flags.
# -Wl,--strip-all:     remove any symbol information from the final binary.

CXXFLAGS_STRIPDEBUG:=

ifeq ($(FLAG_CXX),GCC)
ifeq ($(FLAG_DO_STRIP),yes)
	CXXFLAGS_STRIPDEBUG += -ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--strip-all
endif
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
#     Include the overwrite file again,       #
#     this time the strings will be set       #
#                                             #
###############################################

# Use this file to overwrite the default settings above on a local machine.
# Here, also the strings defined earlier in this file can be overwritten
# if that is what you want.
-include $(projectdir)/OVERWRITE.COMPILE.mk


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
CXXFLAGS += ${CXXFLAGS_STRIPDEBUG}

CXXFLAGS := $(strip $(CXXFLAGS))

CXXFLAGS_EXECUTABLE:=
CXXFLAGS_EXECUTABLE+=$(CXXFLAGS)

# Attempt whole program optimization:
# * fwhole-program : the entire program is visible.
# * fuse-linker-plugin : possible if the gold linker is used (currently buggy).
ifeq ($(FLAG_DO_OPTIMIZE),yes)
	ifeq ($(FLAG_CXX),ICC)
		CXXFLAGS_EXECUTABLE+=-fwhole-program
		# CXXFLAGS_EXECUTABLE+=-fuse-linker-plugin
	else
		ifeq ($(FLAG_CXX),GCC)
			CXXFLAGS_EXECUTABLE+=-fwhole-program
			# CXXFLAGS_EXECUTABLE+=-fuse-linker-plugin
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

CURRENT_COMMIT_ID := $(shell git rev-parse HEAD)
CPPFLAGS += -D GIT_COMMIT_ID=\"$(CURRENT_COMMIT_ID)\"

ifeq ($(FLAG_DISABLE_CHECK_MESHES),yes)
CPPFLAGS += -DDO_NOT_CHECK_MESHES
endif

ifeq ($(FLAG_USE_BACKTRACER),yes)
CPPFLAGS += -DUSE_BACKTRACER
endif

ifeq ($(FLAG_DISABLE_ASSERTIONS),yes)
CPPFLAGS += -DNDEBUG
endif

ifeq ($(FLAG_DISCARD_ASSERT_MESSAGES),yes)
CPPFLAGS += -DDISCARD_ASSERT_MESSAGES
endif

ifeq ($(FLAG_USE_ORIGINAL_ASSERT_MACRO),yes)
CPPFLAGS += -DUSE_ORIGINAL_ASSERT_MACRO
endif

ifeq ($(FLAG_USE_PRIMITIVE_LOGGING),yes)
CPPFLAGS += -DUSE_PRIMITIVE_LOGGING
endif

ifeq ($(FLAG_COLORED_OUTPUT),yes)
CPPFLAGS += -DUSE_COLORED_OUTPUT
endif

ifeq ($(FLAG_DISABLE_STDLIBDEBUG),yes)
else
CPPFLAGS += -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_ASSERTIONS -D_GLIBCXX_SANITIZE_VECTOR
endif

ifeq ($(FLAG_DO_USE_EXTENDED_PRECISION),yes)
CPPFLAGS += -DEXTENDED_PRECISION
endif
ifeq ($(FLAG_DO_USE_SINGLE_PRECISION),yes)
CPPFLAGS += -DSINGLE_PRECISION
endif

ifeq ($(FLAG_DO_OPTIMIZE),yes)
CPPFLAGS += -DELIDE_HOT_FUNCTIONS
endif

ifeq ($(FLAG_DO_OPTIMIZE),yes)
CPPFLAGS += -DELIDE_HOT_FUNCTIONS
endif

# If Windows, then we are probably using MinGW32, and switch to the MinGW32 implementation.
ifeq ($(OS),Windows_NT)
CPPFLAGS += -D __USE_MINGW_ANSI_STDIO=1
endif

# If using Clang, request that transitive includes are disabled.
ifeq ($(FLAG_CXX),CLANG)
CPPFLAGS += -D_LIBCPP_REMOVE_TRANSITIVE_INCLUDES
endif

CPPFLAGS := $(strip $(CPPFLAGS))


###############################################
#                                             #
#         LDLIBS - gather linker flags        #
#                                             #
#         - OpenMP                            #
#         - TCMalloc                          #
#                                             #
###############################################

LDLIBS :=

# If desired, use the minimal version of TCMalloc.
ifeq ($(FLAG_USE_TCMALLOC),yes)
	LDLIBS +=-l:libtcmalloc_minimal.so.4
#	LDLIBS +=-ltcmalloc_minimal
#	LDLIBS +=-ltcmalloc
endif

# On Linux, OpenMP requires linking against the atomic library.
ifeq ($(FLAG_ENABLE_OPENMP),yes)
ifneq ($(OS),Windows_NT)
	LDLIBS +=-latomic
endif
endif

ifeq ($(FLAG_USE_BACKTRACER),yes)
	LDLIBS += -lstdc++_libbacktrace
endif


LDLIBS := $(strip $(LDLIBS))













###############################################
#                                             #
#         Choose the type of linking          #
#         - objectfile                        #
#         - static                            #
#         - dynamic                           #
#         - unspecified                       #
#                                             #
###############################################


LINKINGTYPE:=unspecified

ifeq ($(FLAG_DO_OPTIMIZE),yes)
	LINKINGTYPE:=objectfile
else
	ifeq ($(OS),Windows_NT)
		LINKINGTYPE:=static
	else
		LINKINGTYPE:=dynamic
	endif
endif

# LINKINGTYPE:=static















###############################################
#                                             #
#             Show all parameters             #
#                                             #
###############################################

# Print all the compilation flags set manually or automatically.
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
	$(info FLAG_DO_USE_SINGLE_PRECISION   = $(FLAG_DO_USE_SINGLE_PRECISION) )
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
	$(info )
	@true
	$(info $(PATH))
	# echo $(PATH)
	# echo $$(PATH)
	# echo %PATH%
	$(info MAKE                     = $(MAKE) )
	$(info CXX                      = $(CXX) )
	$(info CXXFLAGS_LANG            = $(CXXFLAGS_LANG) )
	$(info CXXFLAGS_OPTIMIZE        = $(CXXFLAGS_OPTIMIZE) )
	$(info CXXFLAGS_CODEGEN         = $(CXXFLAGS_CODEGEN) )
	$(info CXXFLAGS_WARNINGS        = $(CXXFLAGS_WARNINGS) )
	$(info CXXFLAGS_DEBUG           = $(CXXFLAGS_DEBUG) )
	$(info CXXFLAGS_PROF            = $(CXXFLAGS_PROF) )
	$(info SANITIZERS               = $(SANITIZERS) )
	$(info CXXFLAGS_SANI            = $(CXXFLAGS_SANI) )
	$(info CXXFLAGS_MALLOC          = $(CXXFLAGS_MALLOC) )
	$(info CXXFLAGS_STATICANALYSER  = $(CXXFLAGS_STATICANALYSER) )
	$(info CXXFLAGS_DIAGNOSISFORMAT = $(CXXFLAGS_DIAGNOSISFORMAT) )
	$(info CXXFLAGS_STRIPDEBUG      = $(CXXFLAGS_STRIPDEBUG) )
	$(info CXXFLAGS_EXECUTABLE      = $(CXXFLAGS_EXECUTABLE) )
	$(info CPPFLAGS                 = $(CPPFLAGS) )
	$(info LDFLAGS                  = $(LDFLAGS) )
	$(info LDLIBS                   = $(LDLIBS) )
	$(info LINKINGTYPE              = $(LINKINGTYPE) )
	@true
