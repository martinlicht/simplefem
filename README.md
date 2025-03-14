[![Compilation](https://github.com/martinlicht/simplefem/actions/workflows/main.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/main.yml)
[![Base tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml)
[![FEM tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml)

```
                  [ FFFFF EEEEE EEEEE  CCCC    +       +   ]
                  [ F     E     E     C        +       +   ]
                  [ FFFF  EEEE  EEEE  C     +++++++ +++++++]
                  [ F     E     E     C        +       +   ]
                  [ F     EEEEE EEEEE  CCCC    +       +   ]

                 [ C++ Finite element library based on FEEC ]
                    [www.github.com/martinlicht/simplefem]
```

![Cpp](https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square)

FEEC++ is a work-in-progress C++ library for finite element methods in the spirit of finite element exterior calculus (FEEC). 
Its goal is to enable easy prototyping for fundamental research in numerical PDE with a FEEC point-of-view. 

FEEC++ implements finite element differential forms of arbitrary (uniform) polynomial degree over simplicial meshes,
including Whitney forms and Sullivans forms, together with the relevant algebraic and metric operations. 
It comes with linear algebra implementations and supports uniform refinement and longest edge bisection.

FEEC++ builds and runs on Linux, Windows (Cygwin), and MacOS. 
The only necessary prerequisites are a C++14 compiler (such as GCC and Clang) and GNU Make.

Most important features:

- [x] Written in C++
- [x] Minimal dependencies: C++14 compiler and GNU make
- [x] Compiles and runs on Linux, Windows & Cygwin, or Windows & MinGW-w64

Finite element features: 

- [x] Meshes in dimension 1, 2, and 3
- [x] uniform mesh refinement and longest edge bisection
- [x] Whitney and Sullivan k-forms of any polynomial degree in any dimension
- [x] Exterior derivative, traces, exterior and interior products
- [x] Metric linear operations such as mass operator and Hodge star operator 
- [x] Mass matrices with constant or non-uniform coefficients
 
C++ design guidelines:

- [x] C++14 with optional C++17 enhancements
- [x] Minimal dependencies: C++14 compiler and GNU Make
- [x] Minimal requirements and fast compilation
- [x] Any necessary external library included
- [x] Fail-fast philosophy

Planned finite element features:

- [ ] Different bases of finite element differential forms
- [ ] Polynomial multigrid
- [ ] Curved geometries (surfaces supported already)
- [ ] Finite element spaces with non-uniform polynomial degree
- [ ] Duality-based error estimators and adaptive strategies


This codebase aims for facilitating easy proof-of-concept implementations for finite element methods not found in the standard textbooks, with minimal dependencies and easily portable with regard to OS and machine power. 

This project has some explicit **non-goals**: the project does neither aim for massively distributed-memory parallelism nor for peak high-performance computing. 



