[![Compilation](https://github.com/martinlicht/simplefem/actions/workflows/main.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/main.yml)
[![Base tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests.yml)
[![FEM 2D tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests_comp.yml)
[![FEM 3D tests](https://github.com/martinlicht/simplefem/actions/workflows/unittests_3D.yml/badge.svg)](https://github.com/martinlicht/simplefem/actions/workflows/unittests_3D.yml)

```
                  [ FFFFF EEEEE EEEEE  CCCC    +       +    ]
                  [ F     E     E     C        +       +    ]
                  [ FFFF  EEEE  EEEE  C     +++++++ +++++++ ]
                  [ F     E     E     C        +       +    ]
                  [ F     EEEEE EEEEE  CCCC    +       +    ]

                 [ C++ Finite element library based on FEEC ]
                    [www.github.com/martinlicht/simplefem]
```

![Cpp](https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square)
![Cpp](https://img.shields.io/badge/-C++17-deepskyblue?logo=c%2B%2B&style=flat-square)
![Cpp](https://img.shields.io/badge/-C++20-deepskyblue?logo=c%2B%2B&style=flat-square)
![Cpp](https://img.shields.io/badge/-C++23-deepskyblue?logo=c%2B%2B&style=flat-square)

**Welcome to the FEEC++ Project!**

This repository contains the source code for **FEEC++**, a work-in-progress C++ finite element library designed in the spirit of **finite element exterior calculus** (FEEC). The library aims to facilitate easy and rapid prototyping for fundamental research on numerical methods from the FEEC perspective.

**FEEC++** is intended to be versatile and self-contained: It builds and runs on Linux, Windows (Cygwin and MinGW), and macOS. Its only prerequisites are a C++14 compiler (such as GCC or Clang) and GNU Make.

This project is currently in a pre-release phase, with frequent and substantial changes expected.
That even includes the name of project, which recently changed from `simplefem` to `FEEC++`.

### Features of FEEC++  
FEEC++ implements finite element spaces of arbitrary (uniform) polynomial degree over simplicial meshes, including Whitney forms and Sullivan forms. Additionally, it comes with all necessary linear algebra subroutines (dense matrices, sparse matrices, solvers) and a mesh library that supports uniform refinement and longest edge bisection.

Finite element features include:

- [x] Simplicial meshes in dimensions 1, 2, and 3
- [x] uniform mesh refinement and longest edge bisection
- [x] Whitney and Sullivan k-forms of any polynomial degree in any dimension
- [x] Exterior derivative, traces, exterior and interior products
- [x] Metric linear operations such as mass operator and Hodge star operator
- [x] Mass matrices with constant or non-uniform coefficients
- [x] Krylov subspace methods for sparse matrices and block systems
- [x] Parallelized SSOR preconditioner for CSR matrices
- [x] Operator preconditioning for block systems
- [x] VTK-ready output

Planned finite element features:

- [ ] Polynomial multigrid
- [ ] Finite element spaces with non-uniform polynomial degrees
- [ ] Additional spectrally optimized bases of finite element differential forms
- [ ] Curved geometries (surfaces already supported)
- [ ] Duality-based error estimators and adaptive strategies
- [ ] Discontinuous Galerkin and hybridized higher-order methods

### C++ Design Guidelines  

- **C++14 with optional C++17/C++20/C++23 enhancements**
- **Minimal dependencies:** Requires only a C++14 compiler and GNU Make
- **Self-contained:** All necessary external libraries are included
- **Minimal requirements and fast compilation**
- **Fail-fast philosophy**
- **Comprehensive testing:** Includes numerous tests and examples

Moreover, the code is regularly tested against various compiler warnings, linters, sanitizers, and memory checkers.

### Project Scope and Roadmap

The key objective of this finite element software is to enable proof-of-concept implementations of finite element methods that go beyond what is usually found in standard textbooks.

This software prioritizes minimal dependencies and portability across different operating systems and hardware capabilities. It is supposed to remain operational even on a budget laptop whilst being able to leverage any additional hardware.

At its current stage, the project does not target massively distributed-memory parallelism or peak high-performance computing. These features may be considered once key milestones have been achieved.

### Installation and Usage

Simply clone the repository and run `make` within the source directory. 

The tests will be executable, you just run any `.out` file within the directory `tests`. 

Proper installation will be available as a feature once the project is further down the roadmap.

### Contact

martin.licht@epfl.ch
