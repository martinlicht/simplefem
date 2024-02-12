
# TODO general

The complexity of the project is at a point
where reorganization has become the natural
and also necessary next step.

This requires a rewriting of the unit tests module by
module to ensure that all functionality is being tested
without relying on the visual inspection of the output. 
This is true for the lower levels of the library.
The higher level parts require more quality time to 
ensure that the tests are meaningful. What's more,
neither benchmarking nor production are relevant for 
testing functionality, and thus should be separated.

In the wake of this, the framework of the library
should be progressively extended and polished according
to demands of the unit tests. That touches join error
handling, logging, output colors, assertions...

Finally, on a larger scale, the file structure of the
library as well as the makefile setup should be fixed.
The makefile structure is too complicated and at times
difficult to understand and maintain. It would be better
to improve the structure.

None of the above can be done in a day, so it most likely
requires regular grinding in order to get it done.



# (HIGH) Profiling
---------

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


# (HIGH) text output

Linear operators should only output the basic information but no raw data.
This only affects operators with long data: diagonal, dense, sparse 
All the small parameters should be expressed explicitly

The dense matrix and float vector classes should implement extra output methods 
for the raw data. They parameters for the desired precision, whether indices should be printed,
and for format (row/column and matrix/list).

All operators show return the 'this' pointer in their print methods.

# (HIGH) Traces, locally and globally 

[ ] Determine the orientation of a facet induced by its supersimplex 
[ ] Trace and extension: just the inclusion of the spanning set and its transpose 
[ ] Global trace, incl. orientation and boundary behavior

# (HIGH) AFW-Basis of Sullivan forms

# (HIGH) Floating point exact comparisons ersetzen durch Funktion mit expliziter semantik

# (HIGH) Augmented integration in all numerical tests 

Once the numerical tests have been cleaned up, the right-hand side should always be computed
with (optional) augmented integration. There should be a parameter 'r_plus' to control the added
interpolation quality of the right-hand side. Notably, if 'r_plus == 0', then there should be 
a fallback that avoid repeated computation of the mass matrix.
Similarly, the errors should be computed with augmented integration.

# (HIGH) Clean up unit tests for the numerical examples 

[ ] Don't compute the norms of the solutions and the rhs unless necessary 
[ ] Don't use MINRES unless necessary 
[ ] Remove legacy comments, try to localize the structure.

# (HIGH) Output of solver component 

The solver component outputs should all contain the iteration number if possible.

# (MEDIUM) Averaging for Sullivan and Whitney spaces

For each element, you extract the coefficients of the local basis associated with a subsimplex using a Gram-Matrix. 
This works for Sullivan and Whitney bases alike, no difference. 
You can then average according to some scheme, such as:
- [ ] pick arbitrary
- [x] weight uniformly
- [ ] weight by volume 
Thus you can always average from the larger into the smaller space.

# (MEDIUM) preparing for multigrid

# (MEDIUM) clean up DenseMatrix subsystem 

The following modules look reasonable
- [ ] simple solvers 
- [ ] general solvers (Gauss-Jordan, QR, Cholesky -> inverse )
- [ ] simple scalar functions
- [ ] complicated operations (transpose,determinant,tensorproduct)
- [x] readwrite is never used: retire 

# (MEDIUM) Inverse operators via templates 

Use templates for the inverse operators to implement the 'composed operator' behavior.
Determine the type of solver at compile time depending on the operator class.
This requires a unified solver interface.






# (DONE) Change the include orders 

Go from the most general down to the most specific.
This ensures any overwriting of macros stays local.
Within each grouping, sort alphabetically.

# (DONE) Define and adopt a custom assert macro

There is a function that performs the assert, 
and a macro that delivers the line number and file name
to a function invocation. No further frills.
Use the custom assert macro throughout the project.


# (HIGH) TODO file names 

???Correct the file names in the VTK output of the different mesh and fem test programs ???








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














# (HIGH) Argument names in all header files 
    
The function/method declarations in the header files should provide argument names. 
The names should coincide with the ones in the code but that's not necessary. 

Rationale: this improves readability.







# (HIGH) Question: what are best practices to keep the unit tests up to date with the code?











# (HIGH) Precisions for solvers 

The linear algebra solvers may work with any type of precision, 
such as float or double. Replace the 'magic numbers' in the library
by multiples of the machine epsilon. Generally speaking, 
try to find a good stopping criterion.



# (LOW) Rewrite the unit tests for combinatorics

Generally speaking, the combinatorics unit tests should be more excessive. 
Don't shy away from testing basically everything you could possibly think of as relevant. 
Then move on with the same spirit to 'operators' and the other objects.






# (MEDIUM) Style checker and configuration

Include a style checker such as KWstyle 
and add the necessary configuration files 

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

# (MEDIUM) guarded element access 

All objects that feature element access via brackets,
either blocky brackets or round brackets,
also feature an additional at-method with the same effective behavior. 
The difference is that the at-methods 
always perform bound checks,
which may not the case for the bracket access methods.

- Enforce the effective behavior
- Enforce the bound check policy.


# (MEDIUM) Unit test descriptions

Update the unit test **descriptions** in every module. They seem to be off in many regards.


# (MEDIUM) Revise logging output 

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

# (MEDIUM) Introduce a LOG switch 

Make the logging framework optional by introducing a macro switch 
that enables/disables the logging framework

Then introduce the logging framework throughout the entire code	uniformly.

This requires that the logging interface should be used in the same way
as the entire script for the logging stuff.

# (MEDIUM) Logging class 

Even though advanced logging control would be desirable, 
for the time being it is sufficient if the logging capabilities 
are merely present.

- First layer: semantic wrappers for the cpp streams 
- Second layer: advanced logging classes for the alias streams
- Third layer: primitive MACROS that wrap

In the long run, it would be nice to use a logging class 
that allows for prefixes, git version, date, time, etc.

Setting this up will require some careful thinking 
and refactoring of the entire code. 
A reasonable approach would be a replacement
of cout and cerr throughout the entire code 
by new derivations of the stream class
which facilitate more behavior.

In a first step, this is just two streams 
with the some functionality as cout and cerr.

In a second step, more functionality may be added.

The logging classes that I have seen use macros to emulate 
different log streams, their usage looks like 

    ```LOG << "here is a message";```

alternatively, I would like to skip the shift operator alltogether 
and perhaps replace by a macro to read 

    ```LOG "Here is a message";```

The nice thing is that the log messages get accummulated in the data structure 
and only on destruction of the temporary object the message gets actually written
in the actual logging object. Thus one can impose various 
prefixes and postfixes. 

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
more abundantly throughtou the code. There shouldn't be any magic numbers 
and no 

## (LOW) Basic unit tests 
Not much is to be done here but everything should look fine and reasonable.

## (LOW) Utility unit tests 
Not much is to be done here but everything should look fine and reasonable.

## (LOW) Combinatorics unit tests
Go through all the methods and features in the combinatorics module and write tests for that.
Find ways to test everything independent of the screen output.

## (LOW) Operators unit tests 
Go through all the methods of the Float vectors and write tests for that
Find ways to test the composition operators efficiently.

## (LOW) Dense unit tests 
Go through all the methods of the dense matrix class and write tests for that
For each algorithm, write unit tests. You can test the bottom down version of each algorithm first and then the convenient wrappers for each algorithm.

## (LOW) Sparse unit tests 
Go through all the methods of the sparse matrix class and write tests for that.
Go over the composition operators and check that they do not change the outcome of the computations.

## (LOW) Solver unit tests 
Go through all the solvers and write useful convergence tests for each of those.
You can group them by matrix-type.

## (LOW) Mesh unit tests 
The unit tests are okay but should be rewritten to make everything seamless and consistent.

## (LOW) VTK unit tests 
There is not much to be written here.



# (LOW) interesting meshes

Use the US states map from Randy's source code 
and implement it here. Try to find other triangulations 
too and integrate them as examples. 

# (LOW) Reduce dense matrix module to core functionality 

It suffices to have the core functions for dense linear algebra 
present in this module. Perhaps the matrix I/O should be externalized?
Definitely the dense linear algebra IO should be solely string-based
and not assume specifics about the module.

# (LOW) Operators as non-member functions

Check the classes for member operator functions.
Except for some particular special cases, 
= () [] ->
we can and should turn them into non-member operators.



# (LOW) Gerschgorin circles for Dense and Sparse Matrices
  
implement the gerschgorin row/column circles and the corresponding maximal estimates 




  
  
  







# (LOW) Command line interface

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
    
    





















# TODO UNCLEAR UTILITY


# (UNCLEAR) Implement LU decomposition with different strategies 
  
  The LU decomposition needs to be implemented 
  with different pivoting strategies:
  row, column, or full pivot. 

# (UNCLEAR) Different elementary solvers 
  
  Implement solution algorithms for special matrix types:
  - diagonal solve 
  - left/right triagonal solve 
  - unit left/right/ triagonal solve 
  - averages between left and right solves 
  
## (UNCLEAR) warning command 

add a warning function to the core functionality 
so that you can emit warnings whenever the need arises
instead of calling std::cout.

## (UNCLEAR) Command line 

The handling of command line arguments will be facilitated 
by a set of functions/classes written precisely for that purpose.

This should be a mere extractor class
and be written in the C-conforming subset of C++.

## (UNCLEAR) Fixed-size dynamic array and adoption

Define a template class for a dynamically allocated array
whose size cannot be changed after allocation. 
Copy the std::vector interface but do not provide 
resizing and capacity information.

Use that fixed-size array throughout your code whenever appropiate,
replacing the old std::vector variables with the new ones.
This applies in particular to the linear algebra classes.
    
## (UNCLEAR) implement minimalist file stream wrapper 

    ```
    openinputfile( std::string );
    openoutputfile( std::string );
    ```
    
## (UNCLEAR) Matrix Market

Make Matrix market subdirectory fly.
Add unit tests, then use a solver for a unit test. 

## (UNCLEAR) Global Index Type

Replace any occurence of 'int' by a user-defined type 'Index'.
That type should be large enough for your purposes 
and compatible with the STL standard library.
For example,
    typedef std::size_t Index;

## (UNCLEAR) add complete constructor interfaces 

Apply the rule of six and declare all constructors explicitly
even if merely setting them to default. 

## (UNCLEAR) Add Header files 

All targets should depend on the corresponding header files as well.

## (UNCLEAR) LICENSE File and Copyright notice 

Include a license file into your software.

Include necessary copyright information in the header 
of each file for the entire project.

The projects can be used as example for this:
- dune 
- fetk
- ngsolve 
- ????!!
- LifeV
- vtk

A uniform licence structure should be agreed upon 
before mass reproducing the tests. That being said,
you can also include the license information 
at later stages of the project through the use 
of some simple text manipulation programs.

So for the unit tests, it's more important 
that we have a common structure of the test modules 
ready to go. 

## (UNCLEAR) Separate build and tests targets 

Introduce two different targets, one for building 
the object files and libraries, and the other for 
building the test files 

all: build tests

## (UNCLEAR) Redesign source code organization: library files 

The compilation should place no temporary or output files 
in the source directories. Instead, all output should be 
put into a designated 'build' directory.

The different source directories should specify
the various makefile rules but otherwise not specify anything.
In particular, no cleaning is necessary in those directories.

The makefile in each source directory
puts its output into the common build directory.

There is only one cleaning command for the entire build directory.

## (UNCLEAR) Redesign source code organization: test files 

The test codes are maintained in the source directory
but the programs are put into the same directory.

## (UNCLEAR) Makefile with implicit rules 

The makefile has implicit rules for cpp files
which can greatly simplify the entire make process.
So we may replace the handwritten rules 
by the implicitly defined rules in many cases. 
We merely need to specify the compiler flags.

# (UNCLEAR) Container template 

Flesh out the container template and maybe put it out on code review.


## (UNCLEAR) namespaces for the project

package everything into a namespace.


## (UNCLEAR) Smart Pointers

Should smart pointers be employed throughout the library to make it more robust against user malpractice?


#  (UNCLEAR) VTK OUTPUT 

The general philosophy of the VTK module 
is to treat VTK as an output format alone. 
  
- Perhaps conversion between data set classes
- Conversion between non-VTK related formats
- Direct IO with mesh and algebra classes via specialized routines 

Simple legacy format:
http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
Reading only ASCII







# TODO
  
  
  - zeige die ersten zeilen aller header dateien an
    for datei in ./*/*.hpp; do head -n 2 $datei; done
    include guards umbenennen 
  
  - class for one-dimensional meshes: graphs with no abandoned nodes
  
  
  
  
  
  
  
  
  - hash-tabelle 
  
  

  
  
  
  
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

  - Replace assert(false) by a project-specific macro
  - Magic floating point constant in code 

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

## (DONE) Clean out 'cout' references throughout code 

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

