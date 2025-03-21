[![Compilation](https://github.com/martinlicht/simplefem/actions/workflows/main.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/main.yml)
[![Base tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml)
[![FEM 2D tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml)
[![FEM 3D tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests_3D.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests_3D.yml)

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

Welcome to the FEEC++ project!

This repository contains the source code for FEEC++, a work-in-progress C++ finite element library that adheres to the spirit of finite element exterior calculus (FEEC). The library aims to enable easy and rapid prototyping for fundamental research on numerical methods whilst taking the FEEC point-of-view. 

FEEC++ aims to be versatile and self-contained: it builds and runs on Linux, Windows (Cygwin and MinGW), and MacOS.
Its only necessary prerequisites are a C++14 compiler (such as GCC and Clang) and GNU Make.

This project is currently in a pre-release phase. It is actively built but frequent substantial changes are likely. 
That includes the name of project, which changed from `simplefem` to `FEEC++` recently.

FEEC++ implements finite element spaces of arbitrary (uniform) polynomial degree over simplicial meshes, including Whitney forms and Sullivans forms.
In addition, it comes with all necessary linear algebra subroutines and a mesh library that supports uniform refinement and longest edge bisection.

Finite element features:

- [x] Simplicial meshes in dimension 1, 2, and 3
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



