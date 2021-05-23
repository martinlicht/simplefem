


# TODO general

Introduce a check script which reports common errors in your cpp file.
For example,
	Replace assert(false) by a project-specific macro

# Introduce a LOG switch

Make the logging framework optional by introducing a macro switch 
that enables/disables the logging framework

Then introduce the logging framework throughout the entire code	uniformly



# Precisions for solvers 

The linear algebra solvers may work with any type of precision, such as float or double. Replace the 'magic numbers' in the library by multiples of the machine epsilon. Generally speaking, try to find a good stopping criterion.

# Rewrite the unit tests for combinatorics
Generally speaking, the combinatorics unit tests should be more excessive. Don't shy away from testing basically everything you could possibly think of as relevant. Then move on with the same spirit to 'operators' and the other objects.





# Solver printing data structure 

report_startup
report_finish_success
report_finish_fail
report_restart
report_alert
report_breakdown

iteration_is_printable



# TODO short term

- correct the file names in the VTK output of the different mesh and fem test programs 
- if the file exists, add a suffix number to identify it


# TODO long term and comprehensively


## Matrix Market

Make Matrix market subdirectory fly.
Add unit tests, then use a solver for a unit test. 

## Global Index Type

Replace any occurence of 'int' by a user-defined type 'Index'.
That type should be large enough for your purposes 
and compatible with the STL standard library.
For example,
    typedef std::size_t Index;

## add complete constructor interfaces 

Apply the rule of six and declare all constructors explicitly
even if merely setting them to default. 

## Add Header files 

All targets should depend on the corresponding header files as well.

## LICENSE File

Include a license file into your software.

## README File

Include a readme file into your software.

## Separate build and tests targets 

Introduce two different targets, one for building 
the object files and libraries, and the other for 
building the test files 

all: build tests

## Redesign source code organization: library files 

The compilation should place no temporary or output files 
in the source directories. Instead, all output should be 
put into a designated 'build' directory.

The different source directories should specify
the various makefile rules but otherwise not specify anything.
In particular, no cleaning is necessary in those directories.

The makefile in each source directory
puts its output into the common build directory.

There is only one cleaning command for the entire build directory.

## Redesign source code organization: test files 

The test codes are maintained in the source directory
but the programs are put into the same directory.

## Makefile with implicit rules 

The makefile has implicit rules for cpp files
which can greatly simplify the entire make process.
So we may replace the handwritten rules 
by the implicitly defined rules in many cases. 
We merely need to specify the compiler flags.




## Clean out cout references throughout code 

grep 'cout' ./*/*.?pp 
Conduct a clean out of all direct references 
to the standard preinstalled streams throughout
the entire project so that the above call
should only return very few instances in the main code
and otherwise only stuff in test files.

Moreover, consider replacing all the other stuff
by references to clog instead of cout.


## warning command 

add a warning function to the core functionality 
so that you can emit warnings whenever the need arises
instead of calling std::cout.


## Copyright information in the header 

Include necessary copyright information in the header 
of each file for the entire project.

The projects can be used as example for this:
- dune 
- fetk
- ngsolve 
- ????!!
- vtk

A uniform licence structure should be agreed upon 
before mass reproducing the tests. That being said,
you can also include the license information 
at later stages of the project through the use 
of some simple text manipulation programs.

So for the unit tests, it's more important 
that we have a common structure of the test modules 
ready to go. 

## Logging abstraction 

Encapsulate cout, cerr, and clog within wrapper objects 
that delegate the input to those streams. 
You can then extend the stream wrappers at a later stage 

## Unit test framework

Agree to a common style for the unit tests 

## Container template 

Flesh out the container template and maybe put it out on code review.

## Unit test file organization 

Restructure the unit test directory such that 
each topic receives its own directory.
This should closely mimic the structure 
of the source code 

## MatrixMarket testing 

Use some suits to utilize the matrix market classes 

## Style checker and configuration

Include a style checker such as KWstyle 
and add the necessary configuration files 

## Logging improvement 

Once the logging abstraction has been completed 
enhance the logger with different functionality.


## Reduce dense matrix module to core functionality 

It suffices to have the core functions for dense linear algebra 
present in this module. Perhaps the matrix I/O should be externalized?
Definitely the dense linear algebra IO should be solely string-based
and not assume specifics about the module.


## DONE

### Operators as non-member functions

Check the classes for member operator functions.
Except for some particular special cases, 
= () [] ->
we can and should turn them into non-member operators.














## Speculative issues 


### Rename basic to base 

Basic has the wrong connotation, it makes more sense 
to call it base.


### Smart Pointers

Should smart pointers be employed throughout the library
to make it more robust against user malpractice?


### namespaces for the project

package everything into a namespace.


### extend unit tests 

Make the unit tests for extensive for each class.
Test every single functionality for its accuracy.
  
  

  **** TODO ***
  
  
  - zeige die ersten zeilen aller header dateien an
    for datei in ./*/*.hpp; do head -n 2 $datei; done
    include guards umbenennen 
  
  Namespaces
    - introduce a namespace for FEECPP
    - subnamespaces below that level 
  
  - class for one-dimensional meshes: graphs with no abandoned nodes
  
  
  
  
  
  
  
  
  - hash-tabelle 
  
  
  Zwischenziele: 
  - zweidimensionales mesh ausgeben und visualisieren 
  - scalar mass matrix aufstellen 
  - funktion interpolieren 
  - interpolante ausgeben 
  --->  
  
  - solve poisson problem with natural boundary conditions 
  - output of solution data, error measurement  
    - interpolation of functions to P\Lambda: punkte aussuchen, polynome auswerten  
  - randbedingungen
  - 
  
  
  
  
  
  
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
    - Feature liste welche erwünscht ist: near/far future 
    
  unit test layout:
    logstream that is reference to std::cout 
    TEST_DECL_MODULE( str );
    TEST_DECL_CLASS ( str );
    TEST_DECL_BASIC (     );
    TEST_DECL_TOPIC ( str );
    TEST_ANNOUNCE();
    
    
    















