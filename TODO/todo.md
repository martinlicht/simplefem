
# General remarks

The complexity of the project is at a point where reorganization has become the natural and also necessary next step.

This requires a rewriting of the unit tests module by module to ensure that all functionality is being tested
without relying on the visual inspection of the output. This is particularly relevant for the lower levels of the library.

The higher level parts require more quality time to ensure that the tests are meaningful. What's more, neither benchmarking
nor production are relevant for testing functionality, and thus should be separated.

In the wake of this, the framework of the library should be progressively extended and polished according to demands
of the unit tests. That touches upon error handling, logging, output colors, assertions...

Finally, on a larger scale, the file structure of the library as well as the makefile setup should be fixed. 
The makefile structure is too complicated and at times difficult to understand and maintain, and improvement is due.

None of the above can be done in a day, so it most likely requires regular grinding in order to get it done.




# HOT FIXES / TINY FIXES

- [x] ensure that self-assignment is handled _explicitly_ whenever assignment operators are defined
- [x] all destructors noexcept 
      <https://clang.llvm.org/extra/clang-tidy/checks/performance/noexcept-destructor.html>
- [x] move constructors and assignments noexcept 
      <https://clang.llvm.org/extra/clang-tidy/checks/performance/noexcept-move-constructor.html>
- [x] all destructors virtual or class is final without parents
- [x] destructors do not check (possible moved-from state)
- [x] Implement power method iteration to get the largest eigenvalue
- [x] Does diffinterpol2D/3D still have strange output?
- [?] What is the numerically stable way to compute the determinant?
- [!] The sqrt of a result subject to more than a machine epsilon error
      will be subject to more than the sqrt of the machine epsilon error
- [x] Handle self-assignment in operator=
- [x] `CONSTANT_FLOATINGPOINT_DATATYPE`: uses the same macro case distinction as the other ones
- [x] remove includes from the definition of .cxx files; test the .cxx files
- [x] enable compilation in different C++ modes and remove the error messages.
   * [x] 14
   * [x] 17
   * [x] 20

- [x] convergence tables should handle different precisions, one way or the other:
      best to internally use long double. Requires settling the printf issue
- [x] Clean the printing methods of the convergence tables.
- [x] Check `debug.hpp`
- [x] .cpp file which contains the logging
- [x] fix the vee product; make the tests run in 2D and 3D
- [x] integrate volume forms and scalar fields
- [x] Enable the MinGW printf implementation via: `#define __USE_MINGW_ANSI_STDIO 1`
- [x] std::exp and the like: single precision versions
- [x] correct the computation in the null space test
- [x] tests output into logs
   * [x] create a silent option for each run in the makefile
   * [x] ensure the tests only output into stdout
   * [x] adapt makefile to create a log
- [x] Get it to run on SCITAS
   * [x] login SCITAS shell
   * [x] transfer git repo 
   * [x] Set up the job framework 
- [x] Don't use MINRES whenever you can use another solver
- [x] convergence tables can compute convergence rates
- [x] Don't compute the norms of the solutions and the RHS unless necessary
- [x] Decide what the purpose each test, remove overhead, clarify output strings
   * [x] Sullivan2D
   * [x] Sullivan3D
   * [x] Whitney2D/3D
   * [x] what are poissontransformed old, old2, and the other one? Retire?
   * [x] lshaped? -> these are Maxwell systems. Annotate them. For example, the Lagrange test should reflect simple things and additional overhead
- [x] Clean up the test for the nullspace computation.
- [x] FEM tests check assertions more thoroughly
- [x] Mesh: improve consistency. Include orientation tests in usual tests to save compile time
- [x] VTK: different outputs
- [x] Combinatorics: make tests independent of screen output
- [x] Dense: test everything thoroughly, even up to smaller rounding errors
- [x] Operators: make things independent of screen output
- [x] `mixedsolver.cpp` should test each variant of Hodge-CRM
- [x] Learn about the __SSE__ macro
- [x] How to turn off particular unused variable warnings? -> Blog post
- [x] FEM: *inc -> inc* 
- [x] Coordinate class, rename methods to be get/clone by vertex/dimension
- [x] Check whether the << and >> and bitwise operations are executed on signed integral types
- [!] Compilation error with: no exceptions, optimizations, OpenMP, sanitizers, TCMalloc, stripping, profiling, gold linker: wontfix
- [x] `LOG << "Polynomial degree: " << min_r << " <= " << r << " <= " << max_r << nl;`
- [x] `-Weffc++`: initializer lists and (const) iterators

```c
#define UNUSED_VARIABLE(x) (void)x
__attribute__((unused))
```

- [ ] Warnings in eigenvalue tests to be neutralized 
- [ ] Disabled code should be removed or marked accordingly, same for trivially true conditions
- [ ] Understand GCC options: `-ffold-simple-inlines -fimplicit-constexpr -fno-implement-inlines ? -fvisibility-inlines-hidden` ?












# HIGH: DO THESE NEXT

## (HIGH) Adaptive solution of the Poisson problem (3h)

Combine the primal and mixed formulation for the Poisson Problem in the AFEM test, together with the Hodge star, and compute a posteriori error estimates. 

## (HIGH) Minor FEM rewrite (1h)

- [ ] Summarize files: indexfunctions, polynomialmassmatrix, utilities -> utilities
- [ ] Summarize: global functions

## (HIGH) Algebraic Preconditioners for 3D **READING** (10h)

These preconditioners are intended for the basis blocks such as stiffness and mass matrices. 
Their purpose is speeding up the computation in 3D.
They are possible alternatives to Gauss-Seidel with Eisenstadt and graph coloring.

- [ ] Partitioning into disjoint blocks, then block diagonal preconditioning.
      * Make every node a root node of itself 
      * For each level Loop 
        * For each cluster 
          * Take the one with lowest-rank and extend patch to the ones with who have not been updated before.
      * For each cluster 
        * compute the inverse of that cluster.
        
- [ ] Overlapping blocks, then color the blocks, then Gauss-Seidel. Multiplicative Schwarz / Gauss-Seidel algorithms.
- [ ] Algebraic multigrid. Read the article on the topic.

## (HIGH) AFW-Basis of Sullivan forms **READING**

- [ ] Write about those bases in your article and detail out their construction.
- [ ] Implement these bases.
- [ ] Compare the condition numbers of the bases.
- [ ] Compare the sparsity on standard triangles: rectangular, regular, split-symmetry

## (HIGH) Parallelism for coordinate format (3h)

- [ ] Provide aligned memory allocation via base class
- [ ] Ensure that `FloatVector` and `SparseMatrix` allocate aligned memory 
- [ ] Use aligned memory parallelism for the coordinate format 
























# MEDIUM: DO-ABLE YET LOWER IMPORTANCE

## (MEDIUM) Visualization (4h)

- [ ] Write up instructions on how to use Paraview
- [ ] Include visualization script
- [ ] `lshaped maxwell`: the glyphs have unexplained gaps. Try to fix those. 
- [ ] How stable is the VTK Python interface?

## (MEDIUM) Debug midpoint refinement (10h)




























# UPKEEP: do these every once in a while 

- Warnings about unused parameters and variables 
- Warnings about shadowed parameters and variables 
- Compile with full excessive warnings and then some 
- Compile with full optimization 
- cpplint, clang-tidy

## Interface proof-reading

Pick some component and proofread its interface. 
Assess the function and parameter names, the return types, and the attributes.
If something does not appear right, then make a fix or a TODO note. 

## Rename identifiers 

Follow the guidelines in renaming identifiers to make the code more readable.

## Documentation in the finite element component 

## Floating point exact comparisons are replaced by functions with explicit semantics (2h) **READING** 

Find all magic numbers throughout the code, via clang-tidy. 
Find any instance of `desired_precision` or `machine_epsilson`.
Replace those with semantic tests. 

<https://beta.boost.org/doc/libs/1_68_0/libs/math/doc/html/math_toolkit/float_comparison.html>
Understand the floating-point comparison functions and import them into this project, mutatis mutandis.



































# DENSE 

**READING**

Rearrange and rename everything in the dense matrix module. One suggestion:

- operations: addition, subtraction, scalar multiplication, scalar division, matrix multiplication, transposition, determinant calculation, inverse calculation.
- hard solvers, factorization
- easy solvers
- manipulations (index oriented, not algebra)
- functions to shelve: transpose, skip r/c, det, cofactor, inv, subdet matrix, tensorprod, trace, gerschgorin/eigenvalue, norm

Based on that:

- core matrix class: functions where the ratio computation/output size is small
- manipulations: transpose, deleting rows and columns, tensor product
- simple solvers
- factorizations
- operations: det, cofactor, subdet, inv

Clean up the Dense Matrix component's unit tests.












































# ITERATIVE SOLVERS 

The following requires more quality time but is of no immediate importance. 

## (MEDIUM) Graph coloring and Gauss-Seidel iteration **READING**

- [x] DOF partitioning of CSR Matrices. That is an instance of the graph coloring problem. 
- [ ] Algebraic multigrid
- [ ] Multiplicative Schwarz / Gauss-Seidel algorithms

## (THEORY) Simplify the solver component, CRM in particular **READING**

The project is suffering from too much complexity in the solver component and their unit tests.
There are different variants of the CRM with experimental observation, but theoretical understanding about speed/robustness.
Rather, there should be a template for several variants of the CRM, so that different variations can be studied systematically.
A prerequisite for that endeavor are working notes on the CGM and CRM, even with preconditioners. 
There are different variants for the numerator and denominator of alpha and beta: the latter even has a possible two-term recursion.

- [ ] Simplify the different CR solvers, identify and document the differences, unify via case distinctions.
- [ ] Equip all C++ solvers with preconditioners 
- [ ] Unified solver interface: it makes more sense to set up the solvers as functions instead of classes. 
- [ ] Compare the different Herzog-Soodhalter implementations
- [ ] Verify Herzog-Soodhalter implementation in `systemsparsesolver.cpp` 
- [ ] sparse solvers: implement restart and print modulos, crm with only one abort criterion
- [x] Herzog-Soodhalter is functional, including the sparse variant.
- [x] Fix the whatever solver or retire it.
- [ ] Identify the source for MINRES solver; the other ones are known.
- [ ] Understand the stopping criteria
- [ ] What is the correct variant of CG? (E-Mail Meurant?)
- [ ] Chebyshev iteration.
- [ ] Rewrite the German Wikipedia MINRES method in C-style code.
- [ ] Find another pseudocode for the MINRES method.
- [ ] Create boilerplate FD matrix (square, Laplace, Dirichlet/periodic). Solve that SPD matrix using all available solvers.

## Various

The SSOR preconditioner gives a massive advantage for the CGM and stiffness matrix of the Poisson problem (magnitudes) with a simple choice of SSOR parameter equal 1. Does that carry over to the CRM?

The diagonal preconditioner seems to work well for the CGM and the mass matrix.

The diagonal entries of the Poisson stiffness matrix seem to converge to about 4 as the mesh is refined uniformly. This is compatible with our theoretical scaling estimates. In particular, the diagonal preconditioner will be ineffective here. 

2. Some reading:

- [ ] The solution of Laplacian problems over L-shaped domains with a singular boundary integral method.
- [ ] ON ERROR ESTIMATION IN THE CONJUGATE GRADIENT METHOD AND WHY IT WORKS IN FINITE PRECISION COMPUTATIONS
- [ ] Give another reading to <https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf> 
      Check whether you improve the CGM or CRM (by using only one single loop)
















# LOW: FIGURE OUT HOW TO HANDLE THESE

## (LOW) Avoid compiler warning about allocation sizes 

```cpp
template< typename T, std::enable_if_t<std::is_integral<Integer>::value, bool> = true >
std::size_t ADMISSIBLE_SIZE_CAST( T t )
{
    if( std::is_signed<T>::value ) assert( 0 <= t )
    assert( t < SOME_LIMIT );
    return static_cast<std::size_t>(t);
}
```

## (LOW) Multiplication operator convention to be reconsidered

The sparse matrix classes use the & operator for deep multiplication and do not define the * operator,
which is reserved for shallow operations. Re-assess that convention, possibly change it. 

## (LOW) Operators as non-member functions?

Check the classes for member operator functions.
Except for some particular special cases: `= () [] ->`
We can and should turn them into non-member operators.

## (LOW) Minor improvements of FEM tests (?)

- [x] FEM tests that depend on convergence are not easy to measure, but we can at least test for finiteness
- [ ] Figure out a scalar product implementation that guarantees positive output.
- [ ] Reduce the rounding errors in the mass matrix even further.

## (LOW) Revised semantics for matrix-vector multiplication (2h+Testing)

Specify which operation should be the most basic one and stick with that one.
Probably best to focus on the most simply y = A x instead of the more complex ones. 

## (LOW) Augmented integration in all numerical tests

Once the numerical tests have been cleaned up, the right-hand side should always be computed with (optional) augmented integration.
There should be a parameter 'r_plus' to control the added interpolation quality of the right-hand side.
Notably, if 'r_plus == 0', then there should be a fallback that avoid repeated computation of the mass matrix.
Similarly, the errors should be computed with augmented integration.

## (LOW) Composed Operator mechanism

Recapitulate the interface for the composed operators, which is the reason that the library implements all those pointer methods. 
Is it possible to overload these and contain the infrastructure within those particular files?

## (LOW) Interesting meshes

Use the US states map from Randy's source code and implement it here.
Try to find other triangulations too and integrate them as examples.

## (LOW) Sparse matrix utilities 

- [ ] Enable chaining with the sort_and_compress function, while having rvalue correctness. Generally speaking, all `mutable` methods should have such chaining.
- [ ] SparseMatrix, CSR-Matrix. Statistics substructure that provides the min-average-max number of off-diagonal elements.
- [ ] Check symmetry, measure deviation from symmetry. Assume sort+compress has been applied.
- [ ] Sorting and compressing the sparse matrices does not take that much time. No output annotations needed.
- [ ] Unit tests to check for sparse matrix sorting and compression. Test those sparse matrix routines with random input.

## (LOW) Solver printing data structure

The iterative solvers should be provided a printing data structure
that describes the desired level of printing.
This object can be constructed in various ways.
Whatever the implementation, it provides semantics for telling
what is supposed to be reported.

```cpp
bool report_startup();
bool report_finish_success();
bool report_finish_fail();
bool report_restart();
bool report_alert();
bool report_breakdown();

bool iteration_is_printable();
```

## (LOW) Warm restarts for system solvers

As for the system solver components, the inner iteration seems to work well and is not of concern. In fact, the warm internal restarts are sufficient to keep the iteration number very low. The mass matrix is not a problem.

However, the outer iteration is insufficiently understood. At this point, we cannot rely on any system preconditioners because those may change the matrix structure.

Q: what happens if warm restarts are disabled in the inner iteration?

## (LOW) Preconditioners to implement

Implement the following as classical iterative solvers:

- [ ] Jacobi preconditioner
- [ ] different scaling preconditioners
- [ ] Gauss-Seidel preconditioner
- [x] SOR+SSOR preconditioner
- [ ] block diagonal preconditioner
- [ ] block Gauss-Seidel preconditioner
- [ ] adjustable Gauss-Seidel preconditioner
- [ ] Polynomial preconditioners

## (LOW) Provide Preconditioned variants for all iterative methods

For each iterative method there should be a preconditioned method available.
New iterative methods should only be added if the preconditioned variant is added too.

## (LOW) Rewrite algorithms to be complex number stable

All algorithms should be written in a manner that is also correct when using complex numbers.
This should be accompanied by a written exposition of Krylov subspace methods.












# INACTIVE UNTIL FURTHER NOTICE

## (INACTIVE) Unreal Engine 

- Revisit the Unreal Engine for ideas.
- https://www.gamedev.net/forums/topic/704525-3-quick-ways-to-calculate-the-square-root-in-c/

## (INACTIVE) Convergence of iterative methods 

How to measure the convergence of iterative system solvers? How about the convergence of finite element methods?

## (INACTIVE) Canonicalize the mass matrix on the go

The commutativity is not satisfied sufficiently.
That seems to be due to the mass matrix, since canonicalization reduces that effect.
Can we canonicalize everything already in the matrix assembly?

## (INACTIVE) Inverse operators via templates

Use templates for the inverse operators to implement the 'composed operator' behavior.
Determine the type of solver at compile time depending on the operator class.
This requires a unified solver interface.

## (INACTIVE) Implement LU decomposition with different strategies

Implement LU decomposition with different pivoting strategies: row, column, or full pivot.

## (INACTIVE) Iterative Methods to implement

The following iterative solvers can be implemented.

- [x] Residual Minimizing Descent
- [x] Conjugate Residual Method
- [x] Conjugate Residual Method on Normal Equations
- [ ] Richardson iteration
- [ ] Gradient energy descent
- [ ] Gradient residual descent
- [ ] Symmetric Lanczos minimum residual method
- [ ] GMRES

## (INACTIVE) Global Index Type

Replace any occurrence of 'int' by a user-defined type 'Index'.
That type should be large enough and compatible with the STL standard library.
Possible definition: `typedef std::size_t Index;`

## (INACTIVE) Signal handlers 

SIGINT handler: in case of abnormal abort, emit the name of test/program.

## (INACTIVE) Elaborate Logging class

Even though advanced logging control would be desirable,
for the time being it is sufficient if the logging capabilities are merely present.

- First layer: semantic wrappers for the C++ streams
- Second layer: advanced logging classes for the alias streams
- Third layer: primitive MACROS that wrap

In the long run, it would be nice to use a logging class that allows for prefixes, git version, date, time, etc.

Setting this up will require some careful thinking and refactoring of the entire code.
A reasonable approach would be a replacement of cout and cerr throughout the entire code
by new derivations of the stream class which facilitate more behavior.

In a first step, this is just two streams with the same functionality as cout and cerr.

In a second step, more functionality may be added.

The logging classes that I have seen use macros to emulate
different log streams, their usage looks like

`LOG << "here is a message";`

alternatively, I would like to skip the shift operator all together
and perhaps replace by a macro to read

`LOG "Here is a message";`

The nice thing is that the log messages get accumulated in the data structure
and only on destruction of the temporary object the message gets actually written
in the actual logging object. Thus, one can impose various prefixes and postfixes.

Encapsulate cout, cerr, and clog within wrapper objects
that delegate the input to those streams.
You can then extend the stream wrappers at a later stage

## (INACTIVE) Rewriting the output of each class

- [x] Every (important) class provides a log method:

    ```cpp
    void lg() const { LOG << *this << std::endl; }
    ```

- [x] Every such class implements the shift pattern as well:

    ```cpp
    ostream& operator<<( T t, ostream& os )
    {
      t.print( os );
    }
    ```

- [x] Print does exactly this:

    ```cpp
    virtual std::string text() const override;

    void print( ostream& os ) const
    {
      os << text() << nl;
    }
    ```

- [ ] The `text()` method only emits 'shallow' data, so that no content of vectors and matrices is shown.
    If such content is to be shown, one can use a method such as

    ```cpp
    std::string fulltext() const;
    ```

    Many of your unit tests will need to be rewritten

## (DONT) Some snippet from linear algebra

```cpp
// TODO: Cholesky with Crout pattern, and other possible patterns
// TODO: Break down condition?

// TODO: Cholesky with Pivoting
void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );
void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );

void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t )
{
    if( t == 0 )
        return;
    if( t == 1 )
        QRFactorization( A, Q, R );
    else {
        DenseMatrix Qw(Q), Qv(Q);
        DenseMatrix Rw(R), Rv(R);
        QRFactorizationRepeated( A, Qw, Rw, t-1 );
        QRFactorization( Qw, Qv, Rv );
        Q = Qv;
        R = Rv * Rw;
    }
}

void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q, unsigned int t )
{
    unreachable();

    if( t == 0 )
        return;
    if( t == 1 )
        LQFactorization( A, L, Q );
    else {
        DenseMatrix Qw(Q), Qv(Q);
        DenseMatrix Lw(L), Lv(L);
        LQFactorizationRepeated( A, Lw, Qw, t-1 );
        LQFactorization( Qw, Lv, Qv );
        Q = Qv;
        L = Lw * Lv;
    }
}
```

## (DONT) mesh.simplicial2D.cpp

```cpp
if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) {

    data_edge_firstparent_triangle[ t_e2 ] = counter_triangles + ot;

} else {

    int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
    while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
    current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
    data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edge;

}

/* TODO: triangle parents of opposing vertex */
if( data_vertex_firstparent_triangle[ t_v2 ] == t_old ) {
    int nextparent_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
    data_vertex_firstparent_triangle[ t_v2 ] = counter_edges;
} else {
    int current_tri = data_vertex_firstparent_triangle[ t_v2 ];
    while( data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] != t_old )
    current_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
    data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] = counter_edge;
}
```

## (DONT) Rewrite core float vector class

Write it up in a manner that is close to the STL vector class.
Notably, you cannot inherit from STL classes:
<https://stackoverflow.com/questions/2034916/is-it-okay-to-inherit-implementation-from-stl-containers-rather-than-delegate>
Using raw pointers seems faster than using the STL for the `FloatVector` class

## (DONT) Implement vector slices

A vector slice refers to a part of a vector. 
The slice knows the original vector and some data determine how to access the original members.
The Best approach would be to introduce an abstract class for vectors that captures the interface. 
Then fork off the original class of vectors and the new slice implementation.

## (DONT) Implement lambda-based vectors

For the most general functionality, the get/set methods can then be given in terms of lambdas that produce the required terms/references on the spot.
Note that read-only vectors can be implemented by having the set operation cause an error.
Alternatively, you can introduce a base class 'readable vector' and then derive your general purpose vector from there.

## (DONT) OpenMP parallelization of Float Vector class

Many of the methods in the float vector class are OpenMP parallelizable.
This will be problematic for smaller vectors, and hence not planned yet.

## (DONT) Parallelize matrix transposition

## (DONT) Smart Pointers

Consider adopting smart pointers for the allocation of sparse matrices because you rarely want to deep-copy these matrices.
That being said, doing so would introduce a notable difference to dense matrices, for which that would be too much overhead.














































# INFRASTRUCTURE 


## (INFRASTRUCTURE) MSVC compilation

Make modifications so that the code also compiles on MSVC. 

- [ ] `_HAS_EXCEPTIONS`: The macro is defined in MSVC. Should I use it in addition to the GCC/Clang macros in the debug component?


## (INFRASTRUCTURE) Reorder the compilation makefile (1h)

Compare with the order of setting CXXFLAGS

- Language options (std, ...)
- OpenMP
- Format of diagnostic settings 
- Warning options 
- Static analysis options 
- Whether debugging information is added `-g`
- Profiling instrumentation 
- Code generation options (whether `-fpic -fno-plt`)
- Optimization flags 
- Macro definitions 
- Linker flags 
- Options if TCMalloc is used 
- Whether to strip debug information: `-ffunction-sections -fdata-sections -Wl,--gc-sections -Wl,--strip-all`

```makefile
CXXFLAGS := 
CXXFLAGS += ${CXXFLAGS_LANG}
CXXFLAGS += ${CXXFLAGS_DIAGNOSISFORMAT}
CXXFLAGS += ${CXXFLAGS_WARNINGS}
CXXFLAGS += ${CXXFLAGS_STATICANALYSER}
CXXFLAGS += ${CXXFLAGS_DEBUG}
CXXFLAGS += $(CXXFLAGS_PROF)
CXXFLAGS += $(CXXFLAGS_SANI)
CXXFLAGS += ${CXXFLAGS_MALLOC}
CXXFLAGS += ${CXXFLAGS_OPTIMIZE}
CXXFLAGS += ${CXXFLAGS_CODEGEN}
```

## (INFRASTRUCTURE) namespaces for the project

For wider usage, pack everything into a project namespace. 

## (INFRASTRUCTURE) Rename include guards

Find a format for the include guards that is in line with common practices. 
Write a script that checks whether these practices are upheld. 

How to show the first 2 lines of all .hpp files: `for datei in ./*/*.hpp; do head -n 2 $datei; done`

## (INFRASTRUCTURE) Question: what are best practices keeping the unit tests up to date with the code?

- The unit test file structure should mirror the source code structure. 
- For each component, the tests should focus on the public interface. 
- Keep the tests independent of the code base. 
- Test-driven development: write the test first and then see to implement the source. 
- Code Coverage tools measure what parts of the code are executed during the tests. 
  *Tools include* gcov+lcov, Bullseye Coverage, Codecov, Clang
- CI/CD integration (long-term)



## (INFRASTRUCTURE) LICENSE File and Copyright notice

Include a license file into your software.

Include necessary copyright information in the header of each file for the entire project.

These projects can be used as example for this:

- dune
- fetk
- ngsolve
- ????!!
- LifeV
- vtk

A uniform license structure should be agreed upon before mass reproducing the tests.
That being said, you can also include the license information at later stages of the project
through the use of some simple text manipulation programs.

So for the unit tests, it's more important to have a common structure ready to go.

## (INFRASTRUCTURE) Style checker and configuration

Include a style checker and add the necessary configuration files:

- [ ] KWstyle: active but not widespread. Offers CDash and CTest integration 
- [x] astyle: useful for smaller projects like this 
- [ ] clang-format
- [ ] uncrustify

```sh
astyle --dry-run --mode=c --options=none --ascii --project=path/to/astylerc --recursive "./*.cpp,*.hpp,*.cxx" --exclude="external" --exclude=".legacy" --exclude=".playground"
```

```file .astylerc
# GENERAL BRACE ATTACHMENT STYLE 
# mozilla, kr, otbs,
--style=mozilla
--attach-namespace
--attach-closing-while

# INDENT STYLE 
--indent=spaces=4
--indent-classes
--indent-labels
--indent-switches
--indent-preproc-block
--indent-preproc-cond
--indent-col1-comments
--min-conditional-indent=0
--max-continuation-indent=120
--indent-lambda

# PADDING OPTIONS
--break-blocks
--pad-oper
--pad-comma
--pad-include=none
--unpad-paren
--pad-paren-in
--unpad-brackets 
--pad-brackets-in

--fill-empty-lines
--squeeze-lines=3

# ALIGNMENT OPTIONS
--align-pointer=type
--align-reference=type

# FORMAT OPTIONS
--attach-return-type
--keep-one-line-statements
--convert-tabs
```

## (INFRASTRUCTURE/ARTICLE) Static analyzer

What static analyzers are available by the different compilers?

```sh
-Wanalyzer-too-complex 
-Wanalyzer-double-fclose 
-Wanalyzer-double-free 
-Wanalyzer-exposure-through-output-file 
-Wanalyzer-file-leak 
-Wanalyzer-free-of-non-heap 
-Wanalyzer-malloc-leak 
-Wanalyzer-possible-null-argument 
-Wanalyzer-possible-null-dereference 
-Wanalyzer-null-argument 
-Wanalyzer-null-dereference 
-Wanalyzer-stale-setjmp-buffer 
-Wanalyzer-tainted-array-index 
-Wanalyzer-unsafe-call-within-signal-handler 
-Wanalyzer-use-after-free 
-Wanalyzer-use-of-pointer-in-stale-stack-frame 
```

? `-fanalyzer-transitivity -fanalyzer-verbosity=level # default`

TODO: write an email to the mailing list about the static analyzer.



## (INFRASTRUCTURE) Makefile improvements 

- [ ] Lint and polish the makefiles. 
- [ ] Include links to the manual in the makefiles for quicker reference.
- [ ] Die automatische dependency generation funktioniert noch nicht. Werden alte dependency Angaben erased?
- [ ] Question: best practices in managing makefiles? What dependencies to make explicit?

## (INFRASTRUCTURE) Makefile with implicit rules

The makefile has implicit rules for.cpp files which can greatly simplify the entire make process.
So we may replace the handwritten rules by the implicitly defined rules in many cases.
We merely need to specify the compiler flags.



## (INFRASTRUCTURE) Redesign source code organization: library files

The compilation should place no temporary, build, or output files in the source directories.
In particular, all files built should be put into a designated 'build' directory.
The makefile in each source directory puts its output into the common build directory.

## (INFRASTRUCTURE) General infrastructure and layout of unit tests / Unit test framework

Agree to a common style for the unit tests. The existing unit tests should be streamlined and polished.
Generally speaking, they should be reduced to tests only: benchmarks should be put into a folder of their own;
examples should be a folder of their own as well. Do not shy away from bringing a few tests out of retirement.
Each test should be included in a wrapper function, the main function provided via inclusion from a unified template.
The header file should forward declare the test function and the strings and pass the command line arguments.
The individual tests should only contain a bool-valued function that takes command line arguments, and a few strings that explain the function.

The point is that you should be able to easily change the Unit test header file for all available tests at once, or maybe integrate it with a unit test framework. Finally, the project should transition to catch2 at some point.

Also study the following unit test frameworks:

- <http://unitpp.sourceforge.net/>
- <https://github.com/burner/sweet.hpp/blob/master/options.hpp>
- <https://github.com/burner/sweet.hpp/blob/master/filesystemtest/filesystemtest.cpp>
- <https://github.com/burner/sweet.hpp/blob/master/fector.hpp>
- <https://github.com/ccosmin/tinytest/blob/master/examples/code1.c>
- <https://github.com/greg-white/sTest>

## (INFRASTRUCTURE) Command line interface

The handling of command line options should be localized. This requires the common format for the unit tests. Preferably written, in the C-conforming subset of C++.
Some ideas for command line options:

- [ ] Print help: --help
- [ ] Control output verbosity: --quiet, --verbosity
- [ ] Output, log and error output locations
- [ ] Parameters such as polynomial degree, depending on the program 

## (INFRASTRUCTURE) Stacktrace on Abort

Implement a pretty stacktrace print using C++23 stacktrace facilities. You must figure out how GCC was configured for compilation.

## (INFRASTRUCTURE) Documentation:

Sphinx seems viable <http://www.sphinx-doc.org/en/stable/>, videos via <https://en.wikipedia.org/wiki/Sphinx_(documentation_generator)>.
This software could be used for personal websites.



## (INFRASTRUCTURE) Replace C++ standard library 

**Custom printf implementation** that handles long doubles and offers extensions.

**Custom string and output library** as a thin wrapper around the C functionality. 

**Fixed-size vector** template whose size cannot be changed after allocation. 
Copy the std::vector interface but do not provide resizing and capacity information.

Upvote: <https://stackoverflow.com/questions/10865957/printf-with-stdstring>

## (INFRASTRUCTURE) What are common formats for meshes in FEM software

- `Gmsh`: <http://gmsh.info/doc/texinfo/gmsh.html>
- `Gambit`: <http://web.stanford.edu/class/me469b/handouts/gambit_write.pdf>
- <https://scicomp.stackexchange.com/questions/23882/what-is-a-common-file-data-format-for-a-mesh-for-fem>

Learn more about `Gmsh` and how to utilize it for this project.

## Linting markdown files

```
sed -i 's/[[:space:]]*$//' filename
https://github.com/markdownlint/markdownlint/blob/main/docs/RULES.md
mdl todo.md -r ~MD009,~MD012,~MD013,~MD026,~MD032,~MD034
``` 




























# DONE!

## (DONE/MEDIUM) Review SVG output for 2D meshes and document / make self-documenting    

## (DONE/MEDIUM) Solver output

Within solvers without explicit output enum, distinguish the following cases of the print modulo:

- *positive* print modulo that value and everything else
- * 0* print everything except convergence steps
- *-1* print only start and exit 
- *else* absolute silence

Any automatic value should be 1/20th of the max iteration count.

## (DONE/MEDIUM) Standard max iteration count 

The maximum iteration count should be only the dimension of the system. 

## (DONE/MEDIUM) Solvers study the size of the RHS

Notably, the class-based solvers should not change any member variables during the solution process.

1. Rename the outside argument to `threshold` and use `tolerance` only internally. (requires minimal change)
2. Compute the internal `tolerance` using the RHS norm.
3. Run all tests and see whether it converges. In tests with random RHSs, ensure to normalize the input 

```cpp
// the argument is now called `_threshold`
Float tolerance = 0.;
for( int i = 0; i < N; i++ ) tolerance += b[i]*b[i];
tolerance = maximum( desired_precision, precision * sqrt(tolerance) );
```

## (DONE/MEDIUM) Solvers print whether they have been successful

```cpp
/* HOW DID WE FINISH ? */
recent_deviation = rMAMr;
if( rMAMr > tolerance ) {
    LOG << "PCRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ")\n";
    LOG << "PCRM process has failed. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;
} else {
    LOG << "PCRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ")\n";
    LOG << "PCRM process has succeeded. (" << recent_iteration_count << "/" << max_iteration_count << ") : " << recent_deviation << "/" << tolerance;

}
```

## (DONE/MEDIUM) Solvers should abort if Ar_r or Ad_r is negative ?

## (DONE) Faster assembly of matrices

For higher-dimensional problems, the matrix assembly is speed up by converting all partial matrices to CSR first, and the combining the CSR matrices.
Furthermore, matrix conjugation seems to speed up the assembly, as was shown by measurements. 

## (DONE) clean up dense matrix subsystem

The following modules look reasonable:

- [x] simple scalar functions into the class
- [x] read and write are never used: retire
- [x] factorizations (Gauss-Jordan, QR, Cholesky -> inverse )
- [x] complicated operations (transpose, determinant, tensorproduct)
- [x] simple solvers

## (DONE) Remove dead code

Make dead code alive again or remove it. Search for the following pieces:

```sh
grep 'if(false' ./*/*pp
grep 'if( false' ./*/*pp
```

## (DONE) 'threshold' should be renamed 'tolerance'

## (DONE) GitHub badge C++ >= 17

- <https://img.shields.io/badge/C++-00599C.svg?style=for-the-badge&logo=C++&logoColor=white>
- <https://img.shields.io/badge/-c++-black?logo=c%2B%2B&style=social>
- <https://img.shields.io/badge/C++-14-blue>
- <https://img.shields.io/badge/-C++14-yellow?logo=c%2B%2B>
- <https://img.shields.io/badge/-C++14-cyan?logo=c%2B%2B&style=flat-square>
- <https://img.shields.io/badge/-C++14-skyblue?logo=c%2B%2B&style=flat-square>
- <https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square>

## (DONE) Feature: Unphysical operations

For testing purposes. Given a broken differential form

- (a) put it into canonical form again.
- (b) Alternative, bring it into an equivalent form randomly.

## (DONE) Basic:

[x] Survey of float
[x] Survey/benchmark of factorial/binomial and unit tests
[x] program with leak, hello world
[x] benchmark for memory allocation

## (DONE) Text output

- Get text operations to solvers and mesh class
- combinatorics, vector & operator classes: emphasize text over print
- text: combinatorics, redirect print
- text: operators, dense, sparse ; redirect print
- text: mesh, redirect print
- redirect through entire code

## (DONE) add complete constructor interfaces

Apply the rule of six and declare all constructors explicitly even if merely setting them to default.

## (DONE) Different elementary solvers

Implement solution algorithms for special matrix types:

- diagonal solve
- left/right triangular solve
- unit left/right/ triangular solve
- averages between left and right solves

## (DONE) Precisions for solvers and magic numbers

The linear algebra solvers may work with any type of precision, such as float or double.
Replace the 'magic numbers' in the library by multiples of the machine epsilon.
Generally speaking, try to find a good stopping criterion.

```sh
grep -E '([0-9]+([eE][-+]?[0-9]+))' ./**/*cpp ./*/*/*pp
grep -E '([-+]?\.[0-9]+([eE][-+]?[0-9]+)?)' ./*/*pp ./*/*/*pp
```

## (DONE) Dependencies

- [x] All targets (tests and modules) depend on the header files as well.
- [x] Moreover, every test depends also on the static/dynamic libraries.

## (DONE) How to abort

Termination is done either via `abort()` or via `throw(0)`, depending on whether exceptions are disabled or not.

## (DONE) Separate build and tests targets

Introduce two different targets, one for building the object files and libraries, and the other for building the test files.

`all: build tests`

## (DONE) VTK OUTPUT

The general philosophy of the VTK module is to treat VTK as an output format alone.

Simple legacy format: <http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf>
Reading only ASCII

## (DONE) Implement Hodge star operation

## (DONE) Rewrite composed operators

Reduce the code complexity of the composed operators using pointers, flags, and templated constructors. 
This reduces code complexity. Then basically retire the clone and heir methods.

## (DONE) enable complex coefficients

Enable complex coefficients in a minimalist fashion: compose complex operators from known (real) operators as composed operator.
This should be implemented as a composed operator.

## (DONE) Retire the commutativity tests

The old test files contain commutativity tests. Those should be retired.

## (DONE) Retire unnecessary Whitney tests

Poisson solvers do not need to be tested for Whitney forms.
It suffices to keep a complicated example with mixed boundary conditions.

## (DONE) Clean out legacy alternative tests in the FEM solver files

Introduce a unit test in solverfem for the Darcy-system; the Maxwell is already there.
First, clean out the non-block systems, then the block systems.
Ensure that everything that is deleted has an analogue in the list of solvers.
The following is recommended:

- cpp Mass: CGM
- csr mass: CGM SSOR
- cpp stiff:
- csr stiff: MINRES-CSR or CGM-SSOR
- systems: Herzog-Soodhalter with operator preconditioning

## (DONE) Introduce a custom check script

Introduce a check script which reports common 'errors' in your .cpp file,
that is, stuff you consider important for the design of your code.
For example,

- [ ] Replace assert(false) by a project-specific macro
- [ ] Magic floating point constant in code

```sh
grep --line-number --recursively --color 'assert(' ./*pp
grep --line-number --recursively --color 'cout' ./*pp
grep --line-number --recursively --color '.*[0-9]' ./*pp
grep --line-number --recursively --color '[0-9]e' ./*pp
```

## (DONE) OpenMP pragmas conditional compilation

Every occurrence of 'pragma omp' should be included with a conditional compilation.
This ensures that no compiler warnings about 'unknown pragmas' are issued when you
compile the code with OpenMP disabled.

## (DONE) Copy assignment operator for mesh classes

Since the copy constructor is defined, there should also be an assignment operator
for the different mesh classes, or it should be deleted explicitly.

## (DONE) Clean out 'cout' references throughout code

`grep 'cout' ./*/*.?pp`
Conduct a clean out of all direct references
to the standard preinstalled streams throughout
the entire project so that the above call
should only return very few instances in the main code
and otherwise only stuff in test files.

Moreover, consider replacing all the other stuff
by references to clog instead of cout.

## (DONE) Unit test for condition numbers of single element matrices

- For dimensions 1, 2, and 3
- All form degrees and polynomial degrees up to 6
- Construct the single element meshes, build mass and stiffness matrices
- Compute their inverses (and check their products)
- Perform QR algorithm to find the eigenvalues

Requires: regular triangle and tetrahedra

## (DONE) Conditional compilation when OpenMP

Furthermore, if OpenMP is enabled, then you should compile with an inclusion of thread-safe random number generation.

Generally speaking, you should replace explicit instances of 'rand' and 'srand'
by wrapper functions. This makes it easier to switch to different implementations
throughout whenever that becomes necessary. For example:

- random_integer();
- seed_random_integer();

## (DONE) Output of solver component

The solver component prints should all contain the iteration number if possible.
Each print should start with the iteration number, followed by the message class, and then all other info.

- RESTARTED
- BREAKDOWN
- (NOTICE)
- WARNING
- INTERIM

## (DONE) Argument names in all header files

The function/method declarations in the header files should provide argument names.
The names should coincide with the ones in the code, but that's not necessary.

Rationale: this improves readability.

## (DONE) Change include orders

Go from the most general down to the most specific.
This ensures any overwriting of macros stays local.
Within each grouping, sort alphabetically.

## (DONE) Define and adopt a custom assert macro

There is a function that performs the assertion, and a macro that delivers the line number and file name
to a function invocation. No further frills. Use the custom assert macro throughout the project.

## (DONE) Phantom coordinate in 2D mesh output

Add a phantom coordinate to the output of 2D meshes to plot functions.

## (DONE) Unit test descriptions

Update the unit test **descriptions** in every module. They seem to be off in many cases.

## (DONE) Revise logging output

The logging procedure needs to be reworked.

In particular, switch to an encapsulated approach: all classes should have the ability to produce logs of themselves.
That way, you can isolate the problem in just a few methods throughout the code.

Basically, implement the following methods:

- text:        produces a string presentation (no nl)
- print:       outputs the text() into a given stream (with nl)
- << operator: outputs the text() into a given stream (no nl)
- lg:          outputs the text into the log (with nl).
               This function may take a preamble argument

Revert the current design of logging output: there shouldn't be any automatic newlines. 
Instead, re-introduce the newlines in the tests and deactivate the automatic newline in the logging object.

## (DONE) Introduce a LOG switch

Make the logging framework optional by introducing a macro switch that enables/disables the logging framework.

Then introduce the logging framework throughout the entire code uniformly.

This requires that the logging interface should be used in the same way as the entire script for the logging stuff.

## (DONE) guarded element access

All objects that feature element access via brackets, either blocky brackets or round brackets,
also feature an additional at-method with the same effective behavior. The difference is that the at-methods
always perform bound checks, which may not the case for the bracket access methods.

- Enforce the effective behavior
- Enforce the bound check policy.

## (DONE) Improve the iterator interface of `IndexRange` to allow the full scope

## (DONE) Git ID extraction as macro

## (DONE) Improve the file names used for output

While using the program name for the file output is nice, it is better if you can also add a prefix.
That way you can separate the output of different subtasks more easily.

## (DONE) Dynamic library dependencies

The test programs depend on the object files only if static linking is enabled.
Otherwise, they don't.

## (DONE) Combinatorics generate multiindices tests must actually test something

## (DONE) Dense Matrix rewrite, part 1

- [x] simple solvers remain in one file
- [x] reconcile the norms of dense matrices with the norms of vectors
- [x] matrix tensor product should be part of the functions
- [x] Cholesky, QR, and Gauss-Jordan should be one file


## (DONE) ND meshes

Move this into playgrounds together with the unit tests.

## (DONE) static library dependency

## (DONE) Compilation mode with object files

Ensure that you can compile everything if the executables only take the object file themselves.
This is a third compilation mode in addition to static and dynamic libraries. Ensure it all compiles.

## (DONE) Prevent warnings from the external libraries

## (DONE) Go over the manuals of GCC and Clang, add more possible warnings

Introduce a larger amount of warnings. Only use those that are not enabled by default.
However, turn off the following warnings:

- [x] -Wc++98-compat-local-type-template-args
- [x] -Wreserved-identifier
- [x] -Wold-style-cast
- [x] -Wcovered-switch-default

Following this, go over the list of warnings and re-order everything for the sake of consistency.
Check what needs to be retired.

## (DONE) Finish the printing in the nullspace computation

Enable for all nullspace vectors printing for any polynomial degree.

## (DONE) Printing of higher order & Clean up of unit tests

Most routines only print if r == 1. Generalize that.

- [x] enable higher-order printing wherever convenient, and provide higher-order printing. Agree on polynomial degree.
- [x] writeCellVector data: print barycentric 2-forms
- [x] writeCellScalar data: print barycentric n-forms
- [x] Apply uniform format to Darcy, Maxwell, and curl-curl
- [x] Apply uniform format to lshaped?
- [x] Apply uniform format to Poisson

## (DONE) speed up canonicalize?

- [x] Reduce the number of zero entries in the canonicalization
- [x] Adapt the randomization accordingly

## (DONE) Fix finite difference tests

Check for notanumber and go over the different tests to ensure they test something.

## (DONE) Robustness under shuffling combinatorial data

The code assumes at several points that IndexMaps and sigmas are ordered,
at least when the IndexMaps contains only one element. Conclusion:

- Shuffling the multiindices is fine.
- Shuffling the sigmas is more difficult.

## (DONE) shake the coordinates in tests where there is no explicit functions living on them

- [x] fem
- [x] solvers?
- [x] several finite element tests

## (DONE) Either implement LQ factorization or retire it completely

Implement the LQ factorization and test it

## (DONE) Some warnings to process:

Understand why a compiler might warn about weak vtables and how to avoid that issue.
This concerns IndexMap and MultiIndex in particular. See also:
<https://stackoverflow.com/questions/23746941/what-is-the-meaning-of-clangs-wweak-vtables>

This is not particularly relevant for this project.

## (DONE) Averaging for Sullivan and Whitney spaces

For each element, you extract the coefficients of the local basis associated with a subsimplex using a Gram-Matrix.
This works for Sullivan and Whitney bases alike, no difference.
You can then average according to some scheme, such as:

- [x] arbitrary choice
- [x] weight uniformly
- [x] weight by volume

Thus, you can always average from the non-conforming into the conforming space.

## (DONE) dependencies for object file compilation

Ensure that object files have all their dependencies correctly defined.

## (DONE) fix warnings about printf truncation

```sh
/mesh/mesh.simplicial2D.cpp: In function ‘std::string render_number(double, int)’:
./mesh/mesh.simplicial2D.cpp:3286:42: warning: ‘% *.*f’ directive output between 2 and 2147483958 bytes may exceed minimum required size of 4095 [-Wformat-truncation=]
 3286 |     snprintf( str, str_number_of_chars, "% *.*f", lead+1+tail, tail, num);
      |                                          ^~~~~~
./mesh/mesh.simplicial2D.cpp:3286:41: note: assuming directive output of 3 bytes
 3286 |     snprintf( str, str_number_of_chars, "% *.*f", lead+1+tail, tail, num);
      |                                         ^~~~~~~~

In file included from ./basic/.all.cpp:3:
./basic/basic.cpp: In function ‘std::string timestamp2digitalcode(const timestamp&)’:
./basic/basic.cpp:129:36: warning: ‘%*ju’ directive output may be truncated writing between 10 and 20 bytes into a region of size 11 [-Wformat-truncation=]
  129 |     snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
      |                                    ^~~~
./basic/basic.cpp:129:13: note: ‘snprintf’ output between 11 and 21 bytes into a destination of size 11
  129 |     snprintf( digits, fulllength, "%*ju", numdigits, (uintmax_t)t );
      |     ~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
```

## (DONE) Choose the variant of the determinant. Make manual computation a special case

## (DONE) The finite element matrices should allow the case r==0 for volume forms

## (DONE) unit test for LQ factorization

## (DONE) Fix output of generate multiindices

## (DONE/HIGH) All one vector

Ensure that all-one vectors are not used for constant functions. That only works if r=1. Name: constant_one

## (DONE) Streamline float vector vs dense matrix

Align different Frobenius norms of vectors and dense matrices.
Align the entire `DenseMatrix` and `FloatVector` classes as much as only possible.

> ## (DONE) Nullspace filter
> 
> Streamline the matrix construction in the nullspace discussion:
> Develop a general class of methods to filter out the nullspace, either with single matrices or Hodge systems.
> The nullspace tests then merely use that component internally.
> 
> - [x] scalar examples: simplify the construction and unify 
> - [x] Hodge and mixed: simplify the construction and unify
> - [x] scalar case unified algorithm
> 
> The general method looks like this:
> 
>   1. Grab some random vector with unit mass 
>         OPTIONAL Orthogonalize against the previous nullspace vectors 
>   2. Filter out a nullspace vector 
>   5. Orthogonalize against the previous nullspace vectors 
>   3. Check whether the result is a nullspace vector 
>   4. Check whether the result is too small 


## (DONE/INFRASTRUCTURE) Rename basic to 'base' or 'general' or 'common'

Basic has the wrong connotation, it makes more sense to call it 'base', 'common' or 'general'. Possible names include:

- base 
- common 
- commons
- core 
- general
- global
- shared
- std

Survey a few important projects to get a sense of what name you should use for this one.
That will give you a sense of what you should do.

- Examples: base, common, core, general, std
- MFEM:     general
- Feelpp:   core
- Lifev:    core
- ngsolve:  std
- Fenics:   common
- ??!!!!:   base 
- Hermes:   common
- concepts: ...
- <https://en.wikipedia.org/wiki/List_of_finite_element_software_packages>

## (DONE/INFRASTRUCTURE) remove code from parent directory

rename `base/include.hpp` in the main folder, move to base subfolder, and have all exec.s include it

## (DONE/INFRASTRUCTURE) DECOUPLE FROM C++ STANDARD LIBRARY

[x] Remove std::precision from the test files and benchmarks
[x] remove iostream when possible, seems to reduce compile time, by a few seconds
[x] Remove std::fstream when possible
[x] Remove iomanip in each test module
[x] In each module, reduce includes in header files, use forward declarations
[x] Move STL references to where they are needed, away from header to the.cpp file

## (DONE) Renaming

- [x] better convention for names in files.hpp: snake_case
- [x] summation: clean up the code
- [x] simplify sorthack and retire sorthack test 

















# (ARTICLES)

- [x] local functions in C++
- [ ] marathon of sort algorithms: int <, int mod <, int hash <
- [ ] Why does clog go into stderr?
- [ ] How to enable DLLs? You can use a .def file...


- [ ] remove transitive includes (once pull request done)
- [ ] debug macros in LLVM
- [ ] debug macros in libstd++
- [ ] compilation time marathon
- [ ] proofread the documentation and sort entries, write blog entries

- [ ] `-fshort-enums`: enumerations, their size, their values
- [ ] `-fno-plt`
- [ ] `-fvisibility=hidden`
- [ ] find the largest factorial with only divisions
- [ ] how to compile shared objects on Linux
- [ ] sanitizers

- [ ] compilation time marathon









2D pullback:

```
k=1	max		| PF |             = max                          = max
k=2	max * min       | PF |  max        = max * min / min              = max

3D pullback:
k=1	max 		| PF |             = max 			  = max 
k=2	max * mid       | PF |  max 	   = max * mid / min 		  = max/min * mid 
k=3	max * mid * min | PF |  max * mid  = max * mid * min / min * mid  = max 
```

--- PRODUCT NAME:

    Coffeec.org
    Coffeecpp.org
    FEECpp.org
    FEEC++

  http://asp-software.org/www/misv_resources/business-articles/how-to-name-an-app-a-program-a-company-or-a-service/

--- SOFTWARE:

 * Provide software for rapid prototyping of tensor-valued finite element methods over simplicial grids.
 * Gain practical experience in the implementation of FEEC 
 
 * YES: Maintainability, Flexibility, Produce qualitative results
 * NO: Performance, industry

 * Future features:
 * * mesh refinement in several dimensions (10h)
 * * better preconditioners (4h)
 * * parallelism (lots and lots)
 








```
GCC
    real    2m12.580s
    user    1m57.733s
    sys     0m14.574s
Clang 
    real    2m4.502s
    user    1m52.565s
    sys     0m11.626s
Clang w/o transitive includes
    real    2m4.510s
    user    1m52.688s
    sys     0m11.574s
```



real    3m41.789s
user    3m24.552s
sys     0m17.091s

real    4m43.761s
user    4m21.795s
sys     0m21.579s





# TOPIC: Paper

1. The different bases and spanning sets for spaces on simplices; how to convert between them, embed them, and reduce to them (averaging)
   How to elevate or reduce the polynomial degree
2. Lagrangian interpolation in different bases
3. The exterior derivative in different bases
4. The wedge in different bases
5. The trace in different bases
6. The mass product in different bases
7. The vee product in different bases
8. The Hodge star in different bases




