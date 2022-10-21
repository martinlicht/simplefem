[![C/C++ CI](https://github.com/martinlicht/simplefem/actions/workflows/main.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/main.yml)
[![Unit tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml)

                [ a short description of what you are doing here]

                    [                                    ]
                    [   PUT LOG HERE WITH SOME ASCII     ]
                    [   ART. FIGURE OUT A NICE NAME TOO  ]
                    [                                    ]

                    [www.github.com/martinlicht/simplefem]

FEECPP is a C++ library for finite element methods in the spirit of finite element exterior calculus (FEEC). 
Its goal is to enable easy prototyping for fundamental research in numerical PDE with a FEEC point-of-view. 

FEECPP implements finite element differential forms of arbitrary (uniform) polynomial degree over simplicial meshes,
including Whitney forms and Sullivans forms, together with the relevant algebraic and metric operations. 
It comes with linear algebra implementations and supports uniform refinement and longest edge bisection.

FEECPP builds and runs on Linux, Windows (Cygwin), and MacOS. 
Its only necessary prerequisites are a C++14 compiler (such as GCC and Clang) and GNU Make.

Most important features:

 - [x] Written in C++
 - [x] Minimal dependencies: C++14 compiler and GNU make
 - [x] Compiles and runs on Linux, Windows+Cygwin, or Windows+MinGW-w64
 - [x] Meshes in dimension 1, 2, and 3
 - [x] uniform mesh refinement and longest edge bisection
 - [x] Whitney and Sullivan k-forms of any polynomial degree in any dimension
 - [ ] Algebraic linear operations such as exterior derivatives and traces 
 - [ ] Metric linear operations such as mass operator and Hodge star operator 
 - [x] Mass matrices with constant or non-uniform coefficients  
 - [ ] Duality-based error error estimators and adaptive mesh refinement

Planned features:

 - [ ] Different bases
 - [ ] Polynomial multigrid
 - [ ] Curved geometries and non-trivial coefficients

This software aims for easy prototyping new methods that are outside of the standard textbooks.
This software aims to be easy to hack and to be portable with regard to OS and machine power. 
There are also explicit **non-goals** of this project: 
currently, neither high performance nor massive distributed-memory parallelism are within the scope of the project for the time being.  
This software does not aim at running well-known methods as fast as possible.



