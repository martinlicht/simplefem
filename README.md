
== Profiling ==

Several options are available to generate and assess profiling data. one particularly easy option involves Valgrind. 

The compiler should be invoked with the '-g' option. Then the `callgrind` tool is invoked from command line,

```
valgrind --tool=callgrind [callgrind options] your-program [program options]
```

to generate a file `callgrind.out.[pid]` with the profiling data. Then the profiling data can be presented in the GUI tool `kcachegrind`, as in 

```
kcachegrind callgrind.out.[pid]
```
