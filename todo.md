


# TODO long term and comprehensively

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
    
    
    
  
  ************* General layout: 
  
  - Reduce dense matrix module to the core functions of dense matrices
  - Module for Matrix I/O 
















