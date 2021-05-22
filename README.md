                [ a short description of what you are doing here]

                    [                                    ]
                    [   PUT LOG HERE WITH SOME ASCII     ]
                    [   ART. FIGURE OUT A NICE NAME TOO  ]
                    [                                    ]

                    [Put centered website address here]

FEECPP is a C++ library for finite element methods in the spirit of 
Finite Element Exterior Calculus (FEEC). Its goal is to enable rapid prototyping 
for fundamental research in numerical PDE with a FEEC point-of-view. 

FEECPP provides finite element spaces of polynomial differential forms of 
arbitrary (uniform) polynomial degree over simplicial meshes in arbitrary dimension and topology. 
It implements uniform refinement and longest edge bisection in any dimension.

The target audience for this software are researchers in numerical partial differential equations 
who are well-versed in C++ and who want a customizable finite element software.

FEECPP builds and runs on Linux, Windows (Cygwin), and MacOS. Its only necessary 
prerequisites are a C++17 compiler (such as GCC and Clang) and GNU Make for automatic building.

Lastly, high performance and massive parallelism are explicit **non-goals** of this project 
for the time being. 

Profiling
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
