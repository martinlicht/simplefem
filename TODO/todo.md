
# TODO general

The complexity of the project is at a point
where reorganization has become the natural
and also necessary next step.

This requires a rewriting of the unit tests module by
module to ensure that all functionality is being tested
without relying on the visual inspection of the output. 
This is particularly relevant for the lower levels of the library.

The higher level parts require more quality time to 
ensure that the tests are meaningful. What's more,
neither benchmarking nor production are relevant for 
testing functionality, and thus should be separated.

In the wake of this, the framework of the library
should be progressively extended and polished according
to demands of the unit tests. That touches upon error handling, 
logging, output colors, assertions...

Finally, on a larger scale, the file structure of the
library as well as the makefile setup should be fixed.
The makefile structure is too complicated and at times
difficult to understand and maintain. It would be better
to improve the structure.

None of the above can be done in a day, so it most likely
requires regular grinding in order to get it done.




# (DONE) Go over the manuals of GCC and Clang, add more possible warnings 

Introduce a larger amount of warnings. Only use those that are not enabled by default. 

However Turn off the following warnings:
- [x] -Wc++98-compat-local-type-template-args
- [x] -Wreserved-identifier
- [x] -Wold-style-cast
- [x] -Wcovered-switch-default

Following this, go over the list of warnings and re-order everything for the sake of consistency. 
Check what needs to be retired.

# (HIGH) orientation tests must be included in usual tests

# (HIGH) Solverfem: options

Streamline the main loop in the different solverfem tests to reduce code redundancy

# (HIGH) Implement LQ factorization or retire it completely

Implement the LQ factorization and test it

# (HIGH) Some warnings to process:

Understand why a compiler might warn about weak vtables and how to avoid that issue. 
This concerns IndexMap and MultiIndex in particular. 
https://stackoverflow.com/questions/23746941/what-is-the-meaning-of-clangs-wweak-vtables

# (HIGH) fix warnings about printf truncation 

/mesh/mesh.simplicial2D.cpp: In function ‘std::string render_number(double, int)’:
./mesh/mesh.simplicial2D.cpp:3286:42: warning: ‘% *.*f’ directive output between 2 and 2147483958 bytes may exceed minimum required size of 4095 [-Wformat-truncation=]
 3286 |     snprintf( str, str_number_of_chars, "% *.*f", lead+1+tail, tail, num);
      |                                          ^~~~~~
./mesh/mesh.simplicial2D.cpp:3286:41: note: assuming directive output of 3 bytes
 3286 |     snprintf( str, str_number_of_chars, "% *.*f", lead+1+tail, tail, num);
      |                                         ^~~~~~~~

In file included from ./basic/.all.cpp:3:
./basic/basic.cpp: In function ‘std::string timestamp2digitalcode(const timestamp&)’:
./basic/basic.cpp:129:36: warning: ‘%*ju’ directive output may be truncated writing between 10 and 20 bytes into a region of size 11 [-Wformat-truncation=]
  129 |     snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
      |                                    ^~~~
./basic/basic.cpp:129:13: note: ‘snprintf’ output between 11 and 21 bytes into a destination of size 11
  129 |     snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
      |     ~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# (HIGH) shake the coordinates in tests where there is no explicit functions living on them 
- [ ] meshes 
- [ ] fem 
- [ ] solvers?
- [ ] several finite element tests

# (HIGH) Dense Matrix rewrite, part 2

Finally, rearrange and rename everything in the dense matrix module. One suggestion:
- operations: addition, subtraction, scalar multiplication, scalar division, matrix multiplication, transposition, determinant calculation, inverse calculation.
- hard solvers, factorization
- easy solvers 
- manipulations (index oriented, not algebra)
- functions to shelve: transpose, skip r/c, det, cofactor, inv, subdet matrix, tensorprod, trace, gerschgorin/eigenvalue, norm

Based on that:
- core matrix class: functions where the ratio computation/output size is small 
- manipulations:     transpose, deleting rows and columns, tensor product 
- simple solvers
- factorizations
- operations: det, cofactor, subdet, inv

# (HIGH) Streamline nullspace tests 

Clean up the nullspace computation in the solver component.
Clean up the test for the nullspace computation. 

# (HIGH) FEM rewrite 

- [ ] Summarize: indexfunctions, local.polynomialmassmatrix, utilities -> utilities
- [ ] Summarize: global functions 

# (HIGH) dependencies for object file compilation  

# (HIGH) clean up DenseMatrix subsystem 

The following modules look reasonable
- [ ] simple solvers 
- [ ] general solvers (Gauss-Jordan, QR, Cholesky -> inverse )
- [ ] simple scalar functions
- [ ] complicated operations (transpose,determinant,tensorproduct)
- [x] readwrite is never used: retire 

# (HIGH) Augmented integration in all numerical tests 

Once the numerical tests have been cleaned up, the right-hand side should always be computed with (optional) augmented integration. 
There should be a parameter 'r_plus' to control the added interpolation quality of the right-hand side. 
Notably, if 'r_plus == 0', then there should be a fallback that avoid repeated computation of the mass matrix.
Similarly, the errors should be computed with augmented integration.

# (HIGH) Floating point exact comparisons ersetzen durch Funktion mit expliziter semantik
# (HIGH) Floating-point comparisons

https://beta.boost.org/doc/libs/1_68_0/libs/math/doc/html/math_toolkit/float_comparison.html

Understand the floating-point comparison functions and import them into this project, mutatis mutandis. 

# (HIGH) Fix solvers or fix Wikipedia 

Herzog Soodhalter funktioniert nun, auch die sparse variant.
Vielleicht sollte das umbenannt werden? Scheint mehr Soodhalter zu sein als Herzog 
Whatever solver hingegen kackt ab 

# (HIGH) Clean up unit tests for the numerical examples 

- [ ] Don't compute the norms of the solutions and the rhs unless necessary 
- [ ] Don't use MINRES unless necessary 

# (HIGH) Question: what are best practices to keep the unit tests up to date with the code?

# (HIGH) AFW-Basis of Sullivan forms

# (HIGH) Profiling

Several options are available to generate and assess profiling data. one particularly easy option involves Valgrind. 

The compiler should be invoked with the '-g' option. Then the `callgrind` tool is invoked from command line,

```
valgrind --tool=callgrind [callgrind options] your-program [program options]
```

to generate a file `callgrind.out.[pid]` with the profiling data. Then the profiling data can be presented in the GUI tool `kcachegrind`, as in 

```
kcachegrind callgrind.out.[pid]
```
The target audience for this software are researchers in numerical partial differential equations who are well-versed in C++ and who want a customizable finite element software.

Another alternative is `gprof` as a GUI for profiling data. 

# (HIGH) Floating point exact comparisons ersetzen durch Funktion mit expliziter semantik

# (HIGH) Rename basic to 'base' or 'general' or 'common'

Basic has the wrong connotation, it makes more sense to call it 'base' or 'general'.

Make a survey of a few important projects to get a sense of what name you should use for this one. 
That will give you a sense of what you should do.

Examples: base, common, core, general, std
MFEM: general 
Feelpp: core
Lifev: core 
ngsolve std 
Fenics: common
concepts: ...












# (MEDIUM) Averaging for Sullivan and Whitney spaces

For each element, you extract the coefficients of the local basis associated with a subsimplex using a Gram-Matrix. 
This works for Sullivan and Whitney bases alike, no difference. 
You can then average according to some scheme, such as:
- [ ] pick arbitrary
- [x] weight uniformly
- [ ] weight by volume 
Thus you can always average from the larger into the smaller space.

# (MEDIUM) preparing for multigrid

# (MEDIUM) Style checker and configuration

Include a style checker such as KWstyle or clang-format, and add the necessary configuration files 

Astyle
#astyle --mode=c --options=none --project=path/to/astylerc --ascii --recursive "./*.cpp,*.hpp,*.cxx" --exclude=".playground" --exclude=".legacy"
#--style=mozilla
--attach-namespace 
--indent=spaces=4
--indent-classes
--indent-preproc-cond
--indent-col1-comments
--break-blocks 
--pad-oper 
--pad-comma
--unpad-paren
--pad-paren-in 
--align-pointer=type 
--align-reference=type
--attach-return-type

# (MEDIUM) Solver printing data structure 

The iterative solvers should be provided a printing data structure 
that describes the desired level of printing.
This object can be constructed in various ways. 
Whatever the implementation, it provides semantics for telling 
what is supposed to be reported.

```
bool report_startup();
bool report_finish_success();
bool report_finish_fail();
bool report_restart();
bool report_alert();
bool report_breakdown();

bool iteration_is_printable();
```

# (MEDIUM) Logging class 

Even though advanced logging control would be desirable, 
for the time being it is sufficient if the logging capabilities are merely present.

- First layer: semantic wrappers for the cpp streams 
- Second layer: advanced logging classes for the alias streams
- Third layer: primitive MACROS that wrap

In the long run, it would be nice to use a logging class that allows for prefixes, git version, date, time, etc.

Setting this up will require some careful thinking and refactoring of the entire code. 
A reasonable approach would be a replacement of cout and cerr throughout the entire code 
by new derivations of the stream class which facilitate more behavior.

In a first step, this is just two streams with the some functionality as cout and cerr.

In a second step, more functionality may be added.

The logging classes that I have seen use macros to emulate 
different log streams, their usage looks like 

    ```LOG << "here is a message";```

alternatively, I would like to skip the shift operator alltogether 
and perhaps replace by a macro to read 

    ```LOG "Here is a message";```

The nice thing is that the log messages get accummulated in the data structure 
and only on destruction of the temporary object the message gets actually written
in the actual logging object. Thus one can impose various prefixes and postfixes. 

Encapsulate cout, cerr, and clog within wrapper objects 
that delegate the input to those streams. 
You can then extend the stream wrappers at a later stage 

# (MEDIUM) Unit test framework

Agree to a common style for the unit tests 

The existing unit tests should be streamlined and polished. 
Generally speaking, they should be reduced to tests only:
benchmarks should be put into a folder of their own;
examples should be a folder of their own as well.
Do not shy away from bringing a few tests out of retirement. 

As tests get more complicated, it will pay off to introduce parameters 
more abundantly throughout the code. 











# (LOW) Rewrite the unit tests for combinatorics

Generally speaking, the combinatorics unit tests should be more excessive. 
Don't shy away from testing basically everything you could possibly think of as relevant. 
Then move on with the same spirit to 'operators' and the other objects.

- [ ] Basic 
- [ ] Utility 
- [ ] Combinatorics: make things independent of screen output 
- [ ] Operators: make things independent of screen output 
- [ ] Dense: test everything thoroughly 
- [ ] Sparse: check that composition does not change the outcome 
- [ ] Solver: meaningful convergence tests?
- [ ] Mesh: make everything consistent 
- [ ] VTK

# (LOW) interesting meshes

Use the US states map from Randy's source code 
and implement it here. Try to find other triangulations 
too and integrate them as examples. 

# (LOW) Operators as non-member functions?

Check the classes for member operator functions.
Except for some particular special cases, 
= () [] ->
we can and should turn them into non-member operators.

# (LOW) Command line interface

The handling of command line arguments will be facilitated 
by a set of functions/classes written precisely for that purpose.

This should be a mere extractor class
and be written in the C-conforming subset of C++.

The project comes with unit tests whose behavior can be controlled via commandline.
Generally, there should only be a few commands to describe what is happening.

    --help
    Display a few helpful lines 
    
    --verbose
    output as much as possible 
    
    --quiet 
    Only output warnings or errors 
    
    --logfile
    specify the file were the logging should be directed to
    
    --errfile
    specify the file were the logging should be directed to
    
    --outfile
    specify the file were the output should be directed to


# (LOW) Iterative Methods to implement

The following iterative solvers can be implemented.
  - [x] Residual Minimizing Descent
  - [x] Conjugate Residual Method 
  - [x] Conjugate Residual Method on Normal Equations
  - [ ] Richardson iteration 
  - [ ] Gradient energy descent 
  - [ ] Gradient residual descent 
  - [ ] Symmetric Lanczos minimum residual method 

# (LOW) GMRES with Restart 

Implement the generalized minimal residual method
where the search directions are rebuilt from scratch
after a fixed number of iteration vectors have 
been constructed.

# (LOW) Rewrite algorithms to be complex number stable 

All algorithms should be written in a manner that is also correct when using complex numbers. 
This should be accompanied by a written exposition of Krylov subspace methods.

# (LOW) Preconditioners to implement 

  - [ ] Jacobi preconditioner 
  - [ ] different scaling preconditioners
  - [ ] Gauss-Seidel preconditioner
  - [x] SOR + SSOR preconditioner 
  - [ ] block diagonal preconditioner 
  - [ ] block gauss-seidel preconditioner 
  - [ ] adjustable gauss-seidel preconditioner 
  - [ ] Polynomial preconditioners 

# (LOW) Provide Preconditioned variants for all iterative methods

For each iterative method there should be a preconditioned method available.
New iterative methods should only be added if the preconditioned variant is added too.

# (LOW) LICENSE File and Copyright notice 

Include a license file into your software.

Include necessary copyright information in the header of each file for the entire project.

The projects can be used as example for this:
- dune 
- fetk
- ngsolve 
- ????!!
- LifeV
- vtk

A uniform licence structure should be agreed upon before mass reproducing the tests. 
That being said, you can also include the license information at later stages of the project 
through the use of some simple text manipulation programs.

So for the unit tests, it's more important to have a common structure ready to go. 

# (LOW) Global Index Type

Replace any occurence of 'int' by a user-defined type 'Index'.
That type should be large enough for your purposes 
and compatible with the STL standard library.
For example,
    typedef std::size_t Index;

# (LOW) Smart Pointers

Should smart pointers be employed throughout the library to make it more robust against user malpractice?

    








# INACTIVE UNTIL FURTHER NOTICE

# (INACTIVE) openMP parallelization of Float Vector class

Many of the methods in the float vector class are openMP parallelizable. 
- Constructors
- zero, scale
- NOT random -> perhaps use srand?
- scalarproductwith
- norm, maxnorm, lpnorm
- add vectors 

# (INACTIVE) Rewrite core float vector class 

Write it up in a manner that is close to the STL vector class.
Perhaps even make it a descendant of std::vector<Float> and wrap it only thinly.
https://stackoverflow.com/questions/2034916/is-it-okay-to-inherit-implementation-from-stl-containers-rather-than-delegate
There seem to be complications, so it should be delayed until further notice.
There is rather a speed-up if we replace it by generic C++ memory allocation.
In particular, it does not really mesh with later efforts of parallelization. 
Furthermore, it is better to entirely hide the implementation from the user.

# (INACTIVE) Implement vector slices 

A vector slice refers to a part of a vector.
The slice knows the original vector and 
some data determine how to access the original members.

Best approach would be to introduce an abstract class
for vectors that captures the interface. 
Then fork off the original class of vectors 
and the new slice implementation. 
SEE ALSO Implement lambda-based vectors

# (INACTIVE) Container template 

Flesh out the container template and maybe put it out on code review.

# (INACTIVE) Implement lambda-based vectors 

The get/set methods can then be given in terms 
of lambdas that produce the required terms/references 
on the spot. This gives the most general functionality.

Note that read-only vectors can be implemented 
by having the set operation cause an error.
Alternatively, you can introduce a base class 'readable vector'
and then derive your general purpose vector from there.

# (INACTIVE) Fixed-size dynamic array and adoption

Define a template class for a dynamically allocated array
whose size cannot be changed after allocation. 
Copy the std::vector interface but do not provide 
resizing and capacity information.

Use that fixed-size array throughout your code whenever appropiate,
replacing the old std::vector variables with the new ones.
This applies in particular to the linear algebra classes.
    
# (INACTIVE) Makefile with implicit rules 

The makefile has implicit rules for cpp files
which can greatly simplify the entire make process.
So we may replace the handwritten rules 
by the implicitly defined rules in many cases. 
We merely need to specify the compiler flags.

# (INACTIVE) namespaces for the project

package everything into a namespace.

# (INACTIVE) Inverse operators via templates 

Use templates for the inverse operators to implement the 'composed operator' behavior.
Determine the type of solver at compile time depending on the operator class.
This requires a unified solver interface.

# (INACTIVE) Implement LU decomposition with different strategies 
  
The LU decomposition needs to be implemented with different pivoting strategies:
row, column, or full pivot. 

# (INACTIVE) Redesign source code organization: library files 

The compilation should place no temporary or output files in the source directories. 
Instead, all output should be put into a designated 'build' directory.

The different source directories should specify the various makefile rules but otherwise not specify anything.
In particular, no cleaning is necessary in those directories.

The makefile in each source directory puts its output into the common build directory.

There is only one cleaning command for the entire build directory.










# TODO

zeige die ersten zeilen aller header dateien an
for datei in ./*/*.hpp; do head -n 2 $datei; done
include guards umbenennen 
  
class for one-dimensional meshes: graphs with no abandoned nodes
  
hash-tabelle 
  
https://scicomp.stackexchange.com/questions/23882/what-is-a-common-file-data-format-for-a-mesh-for-fem
Gmsh file format: http://gmsh.info/doc/texinfo/gmsh.html
Gambit:           http://web.stanford.edu/class/me469b/handouts/gambit_write.pdf
  
Question:
  - how to manage unit test for a software library?
  - how to manage makefiles? what dependencies to make explicit?

Documentation:
  Sphinx seems viable 
  http://www.sphinx-doc.org/en/stable/
  Videos at 
  https://en.wikipedia.org/wiki/Sphinx_(documentation_generator)
  Could be used for personal websites 
  
todo-listen in jedem modul:
  - liste der klassen und standard punkte 
    - fertig
    - check
    - unittest 
  - liste der methodenpakete und standard punkte 
    - unittest 
  - Feature liste welche erwuenscht ist: near/far future 
  
unit test layout:
  logstream that is reference to std::cout 
  TEST_DECL_MODULE( str );
  TEST_DECL_CLASS ( str );
  TEST_DECL_BASIC (     );
  TEST_DECL_TOPIC ( str );
  TEST_ANNOUNCE();
  
  









# DONE!

# (DONE) Remove dead code 

Make dead code alive again or remove it. Search for the following pieces:

grep 'if(false' ./*/*pp
grep 'if( false' ./*/*pp

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Github badge C++ >= 17

https://img.shields.io/badge/C++-00599C.svg?style=for-the-badge&logo=C++&logoColor=white
https://img.shields.io/badge/-c++-black?logo=c%2B%2B&style=social
https://img.shields.io/badge/C++-14-blue
https://img.shields.io/badge/-C++14-yellow?logo=c%2B%2B
https://img.shields.io/badge/-C++14-cyan?logo=c%2B%2B&style=flat-square
https://img.shields.io/badge/-C++14-skyblue?logo=c%2B%2B&style=flat-square
https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square

# (DONE) Feature: Unphysical operations 

For testing purposes. Given a broken differential form
(a) put it into canonical form again.
(b) Alternative, bring it into an equivalent form randomly 

# (DONE) Basic:

[x] Survey of float 
[x] Survey/benchmark of factorial/binomial and unit tests 
[x] program with leak, hello world 
[x] benchmark for memory allocation 

# (DONE) Text output 

- Get text operations to solvers and mesh class 
- combinatorics, vector & operator classes: emphasize text over print 
- text: combinatorics, redirect print 
- text: operators, dense, sparse ; redirect print 
- text: mesh, redirect print 
- redirect through entire code 

# (DONE) add complete constructor interfaces 

Apply the rule of six and declare all constructors explicitly even if merely setting them to default. 

# (DONE) Different elementary solvers 
  
Implement solution algorithms for special matrix types:
- diagonal solve 
- left/right triagonal solve 
- unit left/right/ triagonal solve 
- averages between left and right solves 
  
# (DONE) Precisions for solvers and magic numbers 

The linear algebra solvers may work with any type of precision, such as float or double. 
Replace the 'magic numbers' in the library by multiples of the machine epsilon. 
Generally speaking, try to find a good stopping criterion.

grep -E '([0-9]+([eE][-+]?[0-9]+))' ./**/*cpp ./*/*/*pp
grep -E '([-+]?\.[0-9]+([eE][-+]?[0-9]+)?)' ./*/*pp ./*/*/*pp

# (DONE) Dependencies

- [x] All targets (tests and modules) depend on the header files as well. 
- [x] Moreover, every test depends also on the static/dynamic libraries. 

# (DONE) How to abort

Termination is done either via `abort()` or via `throw(0)`, depending on whether exceptions are disabled or not. 

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Separate build and tests targets 

Introduce two different targets, one for building 
the object files and libraries, and the other for 
building the test files 

all: build tests

# (DONE) VTK OUTPUT 

The general philosophy of the VTK module is to treat VTK as an output format alone. 
  
Simple legacy format:
http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
Reading only ASCII

# (DONE) Implement Hodge star operation 

# (DONE) Rewrite composed operators 

Reduce the code complexity of the composed operators using pointers, flags,
and templated constructors. This reduces code complexity. 
Then basically retire the the clone and heir methods.

# (DONE) enable complex coefficients 

Enable complex coefficients in a minimalist fashion:
compose complex operators from known (real) operators as composed operator
This should be implemented as a composed operator.

# (DONE) Retire the communtativity tests 

The old test files contain commutativity tests. Those should be retired.

# (DONE) Retire unnecessary Whitney tests 

Poisson solvers do not need to be tested for Whitney forms. 
It suffices to keep a complicated example with miex boundary conditions.

# (DONE) Clean out legacy alternative tests in the FEM solver files

Introduce a unit test in solverfem for the Darcy-system; the Maxwell is already there.
First, clean out the non-block systems, then the block systems.
Make sure that everything that is deleted has an analogue in the list of solvers.
The following is recommend:
- cpp Mass: CGM
- csr mass: CGM SSOR
- cpp stiff: 
- csr stiff: minres csr or CGM SSOR
- systems: Herzog-Soodhalter mit operator preconditioning

# (DONE) Introduce a custom check script

Introduce a check script which reports common 'errors' in your cpp file,
that is, stuff you consider important for the design of your code.
For example,

- [ ] Replace assert(false) by a project-specific macro
- [ ] Magic floating point constant in code 

```
grep --line-number --recursively --color 'assert(' ./*pp
grep --line-number --recursively --color 'cout' ./*pp
grep --line-number --recursively --color '.*[0-9]' ./*pp
grep --line-number --recursively --color '[0-9]e' ./*pp
```

# (DONE) OpenMP pragmas conditional compilation

Every occurence of 'pragma omp' should be included with a conditional compilation.
This ensures that no compiler warnings about 'unknown pragmas' are issued when you
compile the code with openMP disabled.

# (DONE) Copy assignment operator for mesh classes

Since the copy constructor is defined, there should also be an assignment operator
for the different mesh classes, or it should be deleted explicitly.

# (DONE) Clean out 'cout' references throughout code 

grep 'cout' ./*/*.?pp 
Conduct a clean out of all direct references 
to the standard preinstalled streams throughout
the entire project so that the above call
should only return very few instances in the main code
and otherwise only stuff in test files.

Moreover, consider replacing all the other stuff
by references to clog instead of cout.

# (DONE) Abgleichen der Gitterweiten bei solverfem 

# (DONE) Unit test for condition numbers of single element matrices 

For dimensions 1, 2, and 3
All form degrees and polynomial degrees up to 6
Construct the single element meshes, build mass and stiffness matrices 
compute their inverses (and check their products)
perform QR algorithm to find the eigenvalues

Requires: regular triangle and tetrahedra

# (DONE) Conditional compilation when openMP

Furthermore, if openMP is enabled, then you should compile with an inclusion of thread-safe random number generation. 

Generally speaking, you should replace explicit instances of 'rand' and 'srand' 
by wrapper functions. This makes it easier to switch to different implementations 
throughout whenever that becomes necessary.

For example:
- random_integer();
- seed_random_integer();

# (DONE) Output of solver component 

The solver component prints should all contain the iteration number if possible.
Each print should start with the iterartion number, followed by the message class, and then all other info
RESTARTED
BREAKDOWN
---------
(NOTICE)
WARNING
INTERIM

# (DONE) Argument names in all header files 
    
The function/method declarations in the header files should provide argument names. 
The names should coincide with the ones in the code but that's not necessary. 

Rationale: this improves readability.

# (DONE) Change the include orders 

Go from the most general down to the most specific.
This ensures any overwriting of macros stays local.
Within each grouping, sort alphabetically.

# (DONE) Define and adopt a custom assert macro

There is a function that performs the assert, 
and a macro that delivers the line number and file name
to a function invocation. No further frills.
Use the custom assert macro throughout the project.

# (DONE) Phantom coordinate in 2D mesh output

Add a phantom coordinate coordinate to the output of 2D meshes to plot functions

# (DONE) Unit test descriptions

Update the unit test **descriptions** in every module. They seem to be off in many cases.

# (DONE) Revise logging output 

The logging procedure needs to be reworked.

In particular, switch to an encapsulated approach: all classes should have the ability
to produce logs of themselves. That way, you can isolate the problem 
in just a few methods throughout the code.

Basically, implement the following methods:

- text:        produces a string presentation (no nl)
- print:       outputs the text() into a given stream (with nl)
- << operator: outputs the text() into a given stream (no nl)
- lg:          outputs the text into the log (with nl).
               This function may take a preamble argument

Revert the current design of logging output: there shouldn't be
any automatic newlines. Instead, re-introduce the newlines in the tests
and deactive the automatic newline in the logging object.

# (DONE) Introduce a LOG switch 

Make the logging framework optional by introducing a macro switch 
that enables/disables the logging framework

Then introduce the logging framework throughout the entire code	uniformly.

This requires that the logging interface should be used in the same way
as the entire script for the logging stuff.

# (DONE) guarded element access 

All objects that feature element access via brackets,
either blocky brackets or round brackets,
also feature an additional at-method with the same effective behavior. 
The difference is that the at-methods 
always perform bound checks,
which may not the case for the bracket access methods.

- Enforce the effective behavior
- Enforce the bound check policy.

# (DONE) Remove dead code 

grep 'if(false' ./*/*pp
grep 'if( false' ./*/*pp

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Improve the iterator interface of IndexRange to allow the full scope

# (DONE) Git ID extraction as macro

# (DONE) Improve the file names used for output

While using the program name for the file output is nice, it is better if you can also add an additional prefix.
That way you can separate the output of different subtasks more easily.

# (DONE) Dynamic library dependencies 

The test programs depend on the object files only if static linking is enabled.
Otherwise, they don't.

# (DONE) Output program info 

While using the program name for the file output is nice, it is better if you can also add an additional prefix.
That way you can separate the output of different subtasks more easily.

# (DONE) Combinatorics generate multiindices tests must actually test something

# (DONE) Dense Matrix rewrite, part 1

- [x] simple solvers remain in one file
- [x] reconcile the norms of dense matrices with the norms of vectors 
- [x] matrix tensor product should be part of the functions
- [x] Cholesky, QR, and Gauss-Jordan should be one file 


# (DONE) Reorganize the legacy directory with subdirectories

There already are a few subdirectories. Extend that a little bit.

# (DONE) ND meshes

Move this into playgrounds together with the unit tests.

# (DONE) static library dependency 
# (DONE) Compilation mode with object files

Make that you can compile everything if the executables only take the object file themselves. 
This is a third compilation mode in addition to static and dynamic libraries.
Make sure it all compiles. 

# (DONE) Prevent warnings from the external stb libraries

