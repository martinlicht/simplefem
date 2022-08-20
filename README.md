[![C/C++ CI](https://github.com/martinlicht/simplefem/actions/workflows/main.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/main.yml)
[![Unit tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml)

                [ a short description of what you are doing here]

                    [                                    ]
                    [   PUT LOG HERE WITH SOME ASCII     ]
                    [   ART. FIGURE OUT A NICE NAME TOO  ]
                    [                                    ]

                    [Put centered website address here]

FEECPP is a C++ library for finite element methods in the spirit of 
Finite Element Exterior Calculus (FEEC). Its goal is to enable easy prototyping 
for fundamental research in numerical PDE with a FEEC point-of-view. 

FEECPP implements finite element differential forms of arbitrary (uniform) polynomial degree over simplicial meshes. This includes Whitney forms and Sullivans forms. It supports uniform refinement and longest edge bisection.

FEECPP builds and runs on Linux, Windows (Cygwin), and MacOS. Its only necessary 
prerequisites are a C++17 compiler (such as GCC and Clang) and GNU Make.

Most important features:
[x] Written in C++
[x] Minimal dependencies: C++17 compiler and GNU make
[x] Compiles and runs on Linux, Windows+Cygwin, or Windows+MinGW-w64
[x] Whitney and Sullivan k-forms of any polynomial degree in any dimension
[ ] Operations such as exterior derivatives, Mass matrices, Traces, Hodge-star 
[x] Meshes in dimension 1, 2, and 3, supporting uniform refinement and longest edge bisection.
[ ] Residual error estimators and adaptive mesh refinement

Planned features:
[ ] Different bases
[ ] Curved geometries and non-trivial coefficients


What this is not:

[ ] This software does not aim at running well-known methods as fast as possible. It aims at prototyping new methods that are outside of the standard textbooks.
[ ] 


Lastly, high performance and massive parallelism are explicit **non-goals** of this project 
for the time being. 


