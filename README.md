
                  [ FFFFF EEEEE EEEEE  CCCC    +       +   ]
                  [ F     E     E     C        +       +   ]
                  [ FFFF  EEEE  EEEE  C     +++++++ +++++++]
                  [ F     E     E     C        +       +   ]
                  [ F     EEEEE EEEEE  CCCC    +       +   ]

                 [ C++ Finite element library based on FEEC ]

                  [ FFFFF EEEEE EEEEE  CCCC    +       +   ]
                  [ F     E     E     C        +       +   ]
                  [ FFFF  EEEE  EEEE  C     +++++++ +++++++]
                  [ F     E     E     C        +       +   ]
                  [ F     EEEEE EEEEE  CCCC    +       +   ]

                    [www.github.com/martinlicht/simplefem]

![Cpp](https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square)

FEEC++ is a C++ library for finite element methods in the spirit of finite element exterior calculus (FEEC). 
Its goal is to enable easy prototyping for fundamental research in numerical PDE with a FEEC point-of-view. 

FEEC++ implements finite element differential forms of arbitrary (uniform) polynomial degree over simplicial meshes,
including Whitney forms and Sullivans forms, together with the relevant algebraic and metric operations. 
It comes with linear algebra implementations and supports uniform refinement and longest edge bisection.

FEEC++ builds and runs on Linux, Windows (Cygwin), and MacOS. 
Its only necessary prerequisites are a C++14 compiler (such as GCC and Clang) and GNU Make.

Most important features:

 - [x] Written in C++
 - [x] Minimal dependencies: C++14 compiler and GNU make
 - [x] Compiles and runs on Linux, Windows+Cygwin, or Windows+MinGW-w64

Finite element features: 
 - [x] Meshes in dimension 1, 2, and 3
 - [x] uniform mesh refinement and longest edge bisection
 - [x] Whitney and Sullivan k-forms of any polynomial degree in any dimension
 - [x] Exterior derivative
   [ ] Traces
 - [x] Exterior and interior products
 - [x] Metric linear operations such as mass operator and Hodge star operator 
 - [x] Mass matrices with constant or non-uniform coefficients  
 - [ ] Duality-based error error estimators and adaptive strategies
 
C++ design guidelines:

 - C++14 with optional C++17 enhancements
 - Minimal dependencies: C++14 compiler and GNU Make
 - Fast compile time
 - Minimal requirements and fast compile time 
 - Any external library included
 - Fail-fast philosophy

Planned finite element features:

 - [ ] Different bases
 - [ ] Polynomial multigrid
 - [ ] Curved geometries
 - [ ] Finite element spaces with non-uniform polynomial degree

This software aims for easy prototyping new methods that are outside of the standard textbooks.

This software aims to be easy to hack and to be portable with regard to OS and machine power. 

This project has some explicit **non-goals**:
the scope of this code does not cover massively distributed-memory parallelism or peak high-performance computing. It's purpose is to facilitate easy proof-of-concept implementations for methods that are not found in standard textbooks.



