





	
#	----------- UPKEEP & MINOR FIXES -------------

[!] The sqrt of a result subject to more than a machine epsilon error will be subject to more than the sqrt of the machine epsilon error
[ ] convergence tables should handle different precisions, one way or the other: best to internally use long double. Requires settling the printf issue
[ ] use custom printf implementation for that handles long doubles
[x] CONSTANT_FLOATINGPOINT_DATATYPE: uses the same macro case distinction as the other ones
[ ] std::exp and the like: single precision versions
[ ] what is the numerically stable way to compute the determinant? 

[ ] remove commutativity check in the FEM tests 

[ ] augmented integration for error checks
[ ] augmented integration for rhs? First, understand integration better. 
[ ] interpolation points: what are best interpolation points with smallest lebesgue constant?

[ ] adaptive interpolation 

[ ] Handle self-assignment in operator=

[x] correct the computation in the nullspace test 

[x] diffinterpol2D/3D hat komische ausgaben 

[x] Check debug.hpp for errors 

[x] cpp file which contains the logging

[ ] Why does CLOG go into stderr?

[x] tests output into logs 
  [x] create a silent option for each run in the makefile 
  [x] make sure the tests only output into stdout 
  [x] adapt makefile to create a log  		
  
[x] Zum laufen bringen auf SCITAS 
  [x] login SCITAS shell 
  [x] git repo transferieren 
  [x] job framework aufsetzen 
  
[ ] experiments:
  Die Konvergenztabellen sollen immer aussagen was gerade approximiert wird 
  Der Printmodulo soll immer ein 1/20 der max-iteration betragen 

	


# ------------------ Makefile ---------------------------------------------

[ ] de-clutter the makefiles
  include links to the manual in the makefiles for quicker reference 

[ ] makefile structure 
  [x] common.recipe -> common.compile.mk
  [ ] vtkclean -> outputclean 
  [ ] 

[ ] Die automatische dependency generation funktioniert noch nicht 
  werden alte dependency angaben erased?


  

# ----------- Rewrite text output -------------

[x] Jede (bedeutende) Klasse soll eine Log-Funktion verwenden:
  
  ```
  void lg() const { LOG << *this << std::endl; } 
  ```

[x] Jede Klasse soll das Shift interface implementieren 
  
  ```
  ostream& operator<<( T t, ostream& os )
  {
    t.print( os );
  }
  ```

[x] Print soll genau das tun: 
  
  ```
  virtual std::string text() const override; 
  
  void print( ostream& os ) const
  {
    os << text() << nl;
  }
  ```

[ ] text() method 
    The text() shall only emit 'shallow' data,
  so that no content of vectors and matrices is shown
  If such content is to be shown, one can use a method 
  such as 
  ```
  std::string fulltext() const;
  ```
  Many of your unit tests will need to be rewritten 

[ ] Retire die erweiterten print methoden (etwa ausgabe) 

[ ] composed operators hierarchisch ausgeben: und zwar nur thin mittels text 

[ ] Retire the print() methods throughout the code.




# ------------ sparse matrix component -------------------------

[ ] can you return anon structs?

[ ] SparseMatrix, CSR-Matrix
  Wieviele eintraege sind mindest/durchschnittlich/maximal in einer Reihe um die Diagonale?
  Ausgabe mittels statistics substruktur -> return stat struct 

[ ] check symmetry 
  schreibe eine methode um die deviation von der symmetry zu untersuchen
  du kannst annehmen, dass sort+compress angewendet worden ist 
  use a map datastructure to investigate that

[ ] Die Sortierung und Kompression der SparseMatrix braucht nicht soviel Zeit
  Lasse den Output weg 



# ------------- SOLVERS, upkeep ---------------

[ ] vereinfache die verschiedenen CR solver 
  [ ] finde die unterschiede und markiere sie im Code 
  [ ] mache fallunterscheidungen und vereinheitliche 
  [ ] commit und merge 
  

[ ] alle C++ solvers mit preconditioner 

[ ] solver logging 
  make the logging more reasonable, 
  e.g., give out iteration numbers at the beginning

[ ] solver unit tests:
  Die Konvergenztabellen sollen immer den Namen des Solvers enthalten 
  Der Printmodulo soll immer ein 1/20 der max-iteration betragen 

[ ] im gesamten code durchweg desired_precision anwenden
    das kann in der solver klasse so gesetzt werden 

[ ] Unified solver interface: 
  Es macht wenig Sinn, einen solver als Klasse aufzuziehen. 
  Besser altmodisch mit Parametern, eventuell mit speziellen Default-Werten 
  Idee: Zur Not kann man einen Solver dann einer Klasse verpacken 
        mittels einer Template-Konstruktion?
      
  Funktionen sind:
  - Matrix, initial guess, rhs
  - residual
  - max steps, desired precision
  - print modus 
  - Preconditioner?

# ------------- SOLVERS, important stuff ---------------

[ ] Multiplicative Schwarz / Gauss-Seidel algorithms 

  [ ] DOF partitioning
    Schreibe einen Algorithmus welcher zu gegebener CSR-Matrix 
    die DOF partitioniert. Das ist eine Instanz des graph coloring problem.
    Die max. Zahl der Farben kann man nach oben abschaetzen durch die Anzahl
    der DOF, aber allgemein genuegt es, greedy vorzugehen.
    We can restrict to symmetric matrices 
    
  [ ] Schreibe einen solver, welcher das Faerbung verwendet 
    vielleicht mit einer Gauss-Seidel methode,
    fuer eine multiplicative Schwarz methode.
    
  [ ] Teste das for the mass matrix 

[ ] verstehe die stopping criteria for iterative solvers etwas besser 
[ ] was ist die richtige variante von CG? (E-Mail Meurant?)
  
[ ] Chebyshev iteration:
  Der Solver scheint zu funktionieren. 
  Als naechstes den mittleren Versuch erasen, 
  dann den externen archivieren. 
  
[ ] Verstehe den Chebyshev solver von theoretischer Seite 
  und auch die verschiedenen Arten ihn aufzuschreiben
  


# ---------- COMPOSED OPERATORS -----------------

[ ] Rekapituliere das interface der composed operators,
    und weswegen es diese make-pointer methoden gibt. 
  Lassen sich diese moeglicherweise ueberladen,
  so dass du sie nur einmal definieren brauchst?





# ----------- LOW IMPORTANCE -------------

[x] Die Tests spiegeln nicht unbedingt den source tree wieder 

[ ] rename basic.hpp in the main folder,
	move to basic subfolder,
	and have all exec.s include it

[ ] rename 'basic' module 
	e.g. common, base, ....

[ ] Pseudo unit tests umordnen nach woanders 
	basic/ benchmark helloworld leak logging
	
[ ] In MFEM gibt es die Option entweder eine shared oder static library zu builden.
	That should be used in this project as well.
	
[ ] Generally speaking, go over the unreal engine and harvest good ideas
	[x] unreachable, specialized for virtual functions 

Upvote: https://stackoverflow.com/questions/10865957/printf-with-stdstring






#----------- UNIT TESTING -------------

[ ] Rewrite unit tests as follows:
    - use meaningful names for the tests 
  - introduce unit test names as a special variable in each test
  - once that is done, externalize the tests 
  - create a single header file that includes all the setup / teardown code 

[ ] new unit test: solver 
  Erstelle eine 08/15 FD matrix (square, laplace, Dirichlet) 
  Solve that SPD matrix using all available solvers 
  
[ ] new unit test: solver 
  Erstelle eine 08/15 FD matrix (square, laplace, periodisch) 
  Solve that SPD matrix using all available solvers 

[x] Rewrite unit tests:
    - die namen der tests als variable 
    - nachpruefen dass die namen auch sinn ergeben
  
[ ] SIGINT handler einbauen:
      im falle eines abbruchs geben den namen des tests aus 
    
[ ] Have a look at the following unit test frameworks:

  - http://unitpp.sourceforge.net/
  - https://github.com/burner/sweet.hpp/blob/master/options.hpp
  - https://github.com/burner/sweet.hpp/blob/master/filesystemtest/filesystemtest.cpp
  - https://github.com/burner/sweet.hpp/blob/master/fector.hpp
  - https://github.com/ccosmin/tinytest/blob/master/examples/code1.c
  - https://github.com/greg-white/sTest
  

[ ] verstehe die command line options von catch2

[ ] Wie kann man unit tests verwenden, 
    welche eine exception werfen? Vielleicht im Catch2 forum anfragen.

  ------------------------------------------------------------------------------
  TimeOfWonder
  
  Hello,
  I have been working with some handwritten unit tests for a C++ project and would like to refactor those into using catch2 as an alternative. 
  
  Most of my tests work with `assert` macros or by throwing exceptions. If a test fails, the corresponding program just crashes. In principle, the unit tests look as follows:
  
  ```
  #include "mylib.hpp"
  int main(){
    do_stuff_just_pass_on_success_and_crash_on_fail();
  }
  ```
  
  I am wondering how to convert those into catch2 unit tests in a minimally invasive way, so I don't have to refactor the entire library (at once).
  
  What is a your recommendation?
  ------------------------------------------------------------------------------
	
	



# --------------------- Scott Meyers book --------------

[ ] 7: declare destructors virtual in polymorphic base classes
[ ] 9: never call virtual functions during construction or destruction 
[ ] 10: have assignments return reference to *this
[ ] 11: handle self-assignment in operator=

Better make operators non-member
	
#--------------- Layout & Comment --------------

[ ] Layout of class declarations 
  Which constructors are declared and in which order?
  What is the order of the different methods?
  Which should be explicitly deleted or defaulted?
  
  This concerns the classes in 
  - combinatorics
  - meshes 
  - operators, dense, sparse
  
[ ] Gehe allen Klassen durch und kommentiere das grundlegende Interface:
  - custom constructors 
  - the standard six
  - your personal standard interface 
  - everything else 
  Dabei beobachtest du auch die copy und move semantics,
  und lieferst eventuell implementierung nach.
  Falls Abweichungen auftauchen, fuegst du kommentare ein 
  oder behebst die Abweichungen, insb. Move semantics 


[ ] Kommentiere die Combinatorics-Klassen 
[ ] Kommentiere die Combinatorics-Funktionen 

[ ] Gleiche das Interface von Sparse und CSR Matrix ab wo sinnvoll

[ ] Abgleichen der FloatVector und DenseMatrix class interfaces 
    Vereinfache soweit es geht

[ ] Kommentiere FloatVector und LinearOperator

[ ] Kommentiere die verschiedenen einfachen Operatoren 
  und den flag operator 

[ ] Kommentiere die composed operators soweit es geht,
  und lasse eventuell das komplizierte pointer zeugs vorest aus 



#	---------------- DECOUPLE FROM C++ ------------------------------

[ ] Remove precision from the test files and benchmarks 
  [x] remove iostream when possible 
  [ ] seems to reduce compile time, by a few seconds 
  [!] take care of the ./basic test modules 

[x] Remove fstream 

[ ] In each module, reduce includes in header files to bare minimum, use forward declarations

[ ] Move STL references to where they are needed, away from basic to the cpp file 
  
[ ] replace 'assert' by 'Assert'

[ ] Custom Strings library 
  Custom Output library 

[x] Remove iomanip in each test module 

[ ] Change to makefile to save a few seconds and make it more reasonable:
  as in the test folder, use inclusions 



	
	
	 
	
	


























# TODO general

The complexity of the project is at a point
where reorganization has become the natural
and also necessary next step.

This requires a rewriting of the unit tests module by
module to ensure that all functionality is being tested
without relying on the visual inspection of the output. 
This is particularly relevant for the lower levels of the library.

The higher level parts require more quality time to 
ensure that the tests are meaningful. What's more,
neither benchmarking nor production are relevant for 
testing functionality, and thus should be separated.

In the wake of this, the framework of the library
should be progressively extended and polished according
to demands of the unit tests. That touches upon error handling, 
logging, output colors, assertions...

Finally, on a larger scale, the file structure of the
library as well as the makefile setup should be fixed.
The makefile structure is too complicated and at times
difficult to understand and maintain. It would be better
to improve the structure.

None of the above can be done in a day, so it most likely
requires regular grinding in order to get it done.



# (ARTICLES)  

- [ ] -fshort-enums: enumerations, their size, their values 
- [ ] -fno-plt
- [ ] -fvisibility=hidden
- [ ] find the largest factorial with only divisions 
- [ ] 




# (HIGH) _HAS_EXCEPTIONS
The macro is defined in MSVC. Should I use it?

# (HIGH) Compare different Herzog-Soodhalter methods 

# (HIGH) Semantics for matrix vector multiplication 

Specify which operation should be the most basic one and stick with that one. 

# (HIGH) Floating point exact comparisons ersetzen durch Funktion mit expliziter semantik
# (HIGH) Floating-point comparisons

https://beta.boost.org/doc/libs/1_68_0/libs/math/doc/html/math_toolkit/float_comparison.html
Understand the floating-point comparison functions and import them into this project, mutatis mutandis. 

# (HIGH) AFW-Basis of Sullivan forms

- [ ] Write about those bases in your article and detail out their construction 
- [ ] Implement them 
- [ ] Compare the condition numbers of the bases.
- [ ] Compare the sparsity on standard triangles: rectangular, regular, split-symmetry

# (HIGH) Clean up unit tests 

- [x] Don't use MINRES whenever you can use another solver 
- [x] convergence tables can compute convergence rates 
- [x] Don't compute the norms of the solutions and the rhs unless necessary 
- [x] Decide what the purpose each test, remove overhead, clarify output strings
      - [x] Sullivan2D
      - [x] Sullivan3D
      - [x] Whitney2D/3D
      - [x] what are poissontransformed old, old2, and the other one? Retire?
      - [x] lshaped? -> these are maxwell systems. Annotate them
      For example, the Lagrange test should reflect simple things and additional overhead
- [x] Clean up the test for the nullspace computation. 
- [x] FEM tests check assertions more thoroughly
- [x] ~\writing\project.simplefem :      changes all copied 
- [x] ~\writing\project.simplefem.test : changes all copied 
- [x] Mesh: improve consistency. Include orientation tests in usual tests to save compile time 
- [x] VTK: different outputs 
- [x] Combinatorics: make things independent of screen output 
- [x] Dense: test everything thoroughly, even up to smaller rounding errors 
- [x] Operators: make things independent of screen output 
- [ ] Sparse: check that composition does not change the outcome 
- [ ] Solver: meaningful convergence tests?
- [ ] Clean up the nullspace computation in the solver module.
- [ ] include visualization script 
- [ ] Unit tests must check convergence rates 
- [ ] Streamline the main loop in the different solverfem tests to reduce code redundancy
- [ ] mixedsolver should test each variant of Hodge-CRM

# (HIGH) fem/diffelev3D
# (HIGH) the mass matrix suffers from rounding errors. 
# (HIGH) Try to canonicalize on the go

The commutativity is not satisfied sufficiently. 
That seems to be due to the mass matrix, since canonicalization reduces that effect 
Can we canonicalize everything already in the matrix assembly?

# (HIGH) Code Coverage: background and application










# (HIGH) FEM rewrite 

- [ ] Summarize: indexfunctions, polynomialmassmatrix, utilities -> utilities
- [ ] Summarize: global functions 

# (HIGH) Dense Matrix rewrite, part 2

Finally, rearrange and rename everything in the dense matrix module. One suggestion:
- operations: addition, subtraction, scalar multiplication, scalar division, matrix multiplication, transposition, determinant calculation, inverse calculation.
- hard solvers, factorization
- easy solvers 
- manipulations (index oriented, not algebra)
- functions to shelve: transpose, skip r/c, det, cofactor, inv, subdet matrix, tensorprod, trace, gerschgorin/eigenvalue, norm

Based on that:
- core matrix class: functions where the ratio computation/output size is small 
- manipulations:     transpose, deleting rows and columns, tensor product 
- simple solvers
- factorizations
- operations: det, cofactor, subdet, inv




















# (MEDIUM) Solvers study the size of the RHS

```
    // the argument is now called `_threshold`
    Float threshold = _threshold;
    Float RHS_NORM_SQUARED = 0.;    
    for( int i = 0; i < N; i++ ) RHS_NORM_SQUARED += b[i]*b[i];
    threshold *= sqrt(RHS_NORM_SQUARED);
```


# (MEDIUM) Solvers print whether they have been successful 
```
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


# (MEDIUM) Fix Whatever solvers or fix Wikipedia 

- [x] Herzog Soodhalter funktioniert nun, auch die sparse variant.
- [ ] Identify the source for MINRES solver; the other ones are known
- [x] Fix the whatever solver or retire it 

# (MEDIUM) Interesting meshes

Use the US states map from Randy's source code and implement it here. 
Try to find other triangulations too and integrate them as examples. 

# (MEDIUM) Unit test framework

Agree to a common style for the unit tests. 
The existing unit tests should be streamlined and polished. 
Generally speaking, they should be reduced to tests only:
benchmarks should be put into a folder of their own;
examples should be a folder of their own as well.
Do not shy away from bringing a few tests out of retirement. 

As tests get more complicated, it will pay off to introduce parameters 
more abundantly throughout the code. 




# TODO: Solvers should abort if Ar_r or Ad_r is negative ? 

# TODO: Conjugate Residual method: check // rAr = 0.; 

# TODO: Double check Herzog-Soodhalter implementation in systemsparsesolver 

# TODO: embellish the output of complex operators and other operators with optional flag 

# TODO: mergeelementsinsortedlist should be moved into legacy 

# TODO: make general conceptions about how to handle the output of objects in these modules. 

# TODO: test getdimensionclone/loaddimension for the Coordinates class

Generally, we would like to control the output format by some parameter given to each print function. 
We can assume that the parameters belong to some enum class defined within a class declaration and are specifcally tailored to each class. 
They are a purely optional argument for the print method and may be skipped at convenience.


# (LOW) Operators as non-member functions?

Check the classes for member operator functions.
Except for some particular special cases, 
= () [] ->
we can and should turn them into non-member operators.

# (LOW) Gmsh support

Learn more about Gmsh and how to utilize it for this project. 

# (LOW) Style checker and configuration

Include a style checker such as KWstyle or clang-format, and add the necessary configuration files 

Astyle
#astyle --mode=c --options=none --project=path/to/astylerc --ascii --recursive "./*.cpp,*.hpp,*.cxx" --exclude=".playground" --exclude=".legacy"
#--style=mozilla
--attach-namespace 
--indent=spaces=4
--indent-classes
--indent-preproc-cond
--indent-col1-comments
--break-blocks 
--pad-oper 
--pad-comma
--unpad-paren
--pad-paren-in 
--align-pointer=type 
--align-reference=type
--attach-return-type

# (LOW) Logging class 

Even though advanced logging control would be desirable, 
for the time being it is sufficient if the logging capabilities are merely present.

- First layer: semantic wrappers for the cpp streams 
- Second layer: advanced logging classes for the alias streams
- Third layer: primitive MACROS that wrap

In the long run, it would be nice to use a logging class that allows for prefixes, git version, date, time, etc.

Setting this up will require some careful thinking and refactoring of the entire code. 
A reasonable approach would be a replacement of cout and cerr throughout the entire code 
by new derivations of the stream class which facilitate more behavior.

In a first step, this is just two streams with the some functionality as cout and cerr.

In a second step, more functionality may be added.

The logging classes that I have seen use macros to emulate 
different log streams, their usage looks like 

    ```LOG << "here is a message";```

alternatively, I would like to skip the shift operator alltogether 
and perhaps replace by a macro to read 

    ```LOG "Here is a message";```

The nice thing is that the log messages get accummulated in the data structure 
and only on destruction of the temporary object the message gets actually written
in the actual logging object. Thus one can impose various prefixes and postfixes. 

Encapsulate cout, cerr, and clog within wrapper objects 
that delegate the input to those streams. 
You can then extend the stream wrappers at a later stage 

# (LOW) Solver printing data structure 

The iterative solvers should be provided a printing data structure 
that describes the desired level of printing.
This object can be constructed in various ways. 
Whatever the implementation, it provides semantics for telling 
what is supposed to be reported.

```
bool report_startup();
bool report_finish_success();
bool report_finish_fail();
bool report_restart();
bool report_alert();
bool report_breakdown();

bool iteration_is_printable();
```

# (LOW) Command line interface

The handling of command line arguments will be facilitated 
by a set of functions/classes written precisely for that purpose.

This should be a mere extractor class
and be written in the C-conforming subset of C++.

The project comes with unit tests whose behavior can be controlled via commandline.
Generally, there should only be a few commands to describe what is happening.

    --help
    Display a few helpful lines 
    
    --verbose
    output as much as possible 
    
    --quiet 
    Only output warnings or errors 
    
    --logfile
    specify the file were the logging should be directed to
    
    --errfile
    specify the file were the logging should be directed to
    
    --outfile
    specify the file were the output should be directed to


# (LOW) Iterative Methods to implement

The following iterative solvers can be implemented.
  - [x] Residual Minimizing Descent
  - [x] Conjugate Residual Method 
  - [x] Conjugate Residual Method on Normal Equations
  - [ ] Richardson iteration 
  - [ ] Gradient energy descent 
  - [ ] Gradient residual descent 
  - [ ] Symmetric Lanczos minimum residual method 

# (LOW) GMRES with Restart 

Implement the generalized minimal residual method
where the search directions are rebuilt from scratch
after a fixed number of iteration vectors have 
been constructed.

# (LOW) Rewrite algorithms to be complex number stable 

All algorithms should be written in a manner that is also correct when using complex numbers. 
This should be accompanied by a written exposition of Krylov subspace methods.

# (LOW) Preconditioners to implement 

  - [ ] Jacobi preconditioner 
  - [ ] different scaling preconditioners
  - [ ] Gauss-Seidel preconditioner
  - [x] SOR + SSOR preconditioner 
  - [ ] block diagonal preconditioner 
  - [ ] block gauss-seidel preconditioner 
  - [ ] adjustable gauss-seidel preconditioner 
  - [ ] Polynomial preconditioners 

# (LOW) Provide Preconditioned variants for all iterative methods

For each iterative method there should be a preconditioned method available.
New iterative methods should only be added if the preconditioned variant is added too.

# (LOW) LICENSE File and Copyright notice 

Include a license file into your software.

Include necessary copyright information in the header of each file for the entire project.

The projects can be used as example for this:
- dune 
- fetk
- ngsolve 
- ????!!
- LifeV
- vtk

A uniform licence structure should be agreed upon before mass reproducing the tests. 
That being said, you can also include the license information at later stages of the project 
through the use of some simple text manipulation programs.

So for the unit tests, it's more important to have a common structure ready to go. 

# (LOW) Global Index Type

Replace any occurence of 'int' by a user-defined type 'Index'.
That type should be large enough for your purposes 
and compatible with the STL standard library.
For example,
    typedef std::size_t Index;

# (LOW) Smart Pointers

Should smart pointers be employed throughout the library to make it more robust against user malpractice?

    








# INACTIVE UNTIL FURTHER NOTICE

# (INACTIVE): parallelize matrix transposition 

# (INACTIVE) Some snippet from linear algebra 

// TODO: Cholesky with Crout pattern, and other possible patterns 
// TODO: Break down condition?

// TODO: Cholesky with Pivoting 


 // QR repeated 
 // LQ repeated 

// void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );
// void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t );


// // // void QRFactorizationRepeated( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R, unsigned int t )
// // // {
// // //     if( t == 0 )
// // //         return;
// // //     if( t == 1 )
// // //         QRFactorization( A, Q, R );
// // //     else {
// // //         DenseMatrix Qw(Q), Qv(Q);
// // //         DenseMatrix Rw(R), Rv(R);
// // //         QRFactorizationRepeated( A, Qw, Rw, t-1 );
// // //         QRFactorization( Qw, Qv, Rv );
// // //         Q = Qv;
// // //         R = Rv * Rw;
// // //     }
// // // }
// // // 
// // // 
// // // void LQFactorizationRepeated( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q, unsigned int t )
// // // {
// // //     unreachable();
// // //     
// // //     if( t == 0 )
// // //         return;
// // //     if( t == 1 )
// // //         LQFactorization( A, L, Q );
// // //     else {
// // //         DenseMatrix Qw(Q), Qv(Q);
// // //         DenseMatrix Lw(L), Lv(L);
// // //         LQFactorizationRepeated( A, Lw, Qw, t-1 );
// // //         LQFactorization( Qw, Lv, Qv );
// // //         Q = Qv;
// // //         L = Lw * Lv;
// // //     }
// // // }




# (INACTIVE) mesh.simplicial2D.cpp 

//         if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) {
//           
//           data_edge_firstparent_triangle[ t_e2 ] = counter_triangles + ot;
//           
//         } else {
//           
//           int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
//           while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
//             current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
//           data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edge;
//           
//         }
//       
//         /* TODO: triangle parents of opposing vertex */
//         if( data_vertex_firstparent_triangle[ t_v2 ] == t_old ) { 
//           int nextparent_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
//           data_vertex_firstparent_triangle[ t_v2 ] = counter_edges;
//         } else {
//           int current_tri = data_vertex_firstparent_triangle[ t_v2 ];
//           while( data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] != t_old )
//             current_tri = data_triangle_nextparents_of_vertices[ current_tri ][ 2 ];
//           data_triangle_nextparents_of_vertices[ current_tri ][ 2 ] = counter_edge;
//         }


# (INACTIVE) Question: what are best practices to keep the unit tests up to date with the code?

# (INACTIVE) openMP parallelization of Float Vector class

Many of the methods in the float vector class are openMP parallelizable. 
- Constructors
- zero, scale
- NOT random -> perhaps use srand?
- scalarproductwith
- norm, maxnorm, lpnorm
- add vectors 

# (INACTIVE) Rewrite core float vector class 

Write it up in a manner that is close to the STL vector class.
Perhaps even make it a descendant of std::vector<Float> and wrap it only thinly.
https://stackoverflow.com/questions/2034916/is-it-okay-to-inherit-implementation-from-stl-containers-rather-than-delegate
There seem to be complications, so it should be delayed until further notice.
There is rather a speed-up if we replace it by generic C++ memory allocation.
In particular, it does not really mesh with later efforts of parallelization. 
Furthermore, it is better to entirely hide the implementation from the user.

# (INACTIVE) Implement vector slices 

A vector slice refers to a part of a vector.
The slice knows the original vector and 
some data determine how to access the original members.

Best approach would be to introduce an abstract class
for vectors that captures the interface. 
Then fork off the original class of vectors 
and the new slice implementation. 
SEE ALSO Implement lambda-based vectors

# (INACTIVE) Container template 

Flesh out the container template and maybe put it out on code review.

# (INACTIVE) Implement lambda-based vectors 

The get/set methods can then be given in terms 
of lambdas that produce the required terms/references 
on the spot. This gives the most general functionality.

Note that read-only vectors can be implemented 
by having the set operation cause an error.
Alternatively, you can introduce a base class 'readable vector'
and then derive your general purpose vector from there.

# (INACTIVE) Fixed-size dynamic array and adoption

Define a template class for a dynamically allocated array
whose size cannot be changed after allocation. 
Copy the std::vector interface but do not provide 
resizing and capacity information.

Use that fixed-size array throughout your code whenever appropiate,
replacing the old std::vector variables with the new ones.
This applies in particular to the linear algebra classes.
    
# (INACTIVE) Makefile with implicit rules 

The makefile has implicit rules for cpp files
which can greatly simplify the entire make process.
So we may replace the handwritten rules 
by the implicitly defined rules in many cases. 
We merely need to specify the compiler flags.

# (INACTIVE) namespaces for the project

package everything into a namespace.

# (INACTIVE) Inverse operators via templates 

Use templates for the inverse operators to implement the 'composed operator' behavior.
Determine the type of solver at compile time depending on the operator class.
This requires a unified solver interface.

# (INACTIVE) Implement LU decomposition with different strategies 
  
The LU decomposition needs to be implemented with different pivoting strategies:
row, column, or full pivot. 

# (INACTIVE) Redesign source code organization: library files 

The compilation should place no temporary or output files in the source directories. 
Instead, all output should be put into a designated 'build' directory.

The different source directories should specify the various makefile rules but otherwise not specify anything.
In particular, no cleaning is necessary in those directories.

The makefile in each source directory puts its output into the common build directory.

There is only one cleaning command for the entire build directory.

# (INACTIVE) fem/diffinterpol

Can these errors be reduced? Probably not, it's a standard procedure.

# (INACTIVE) Augmented integration in all numerical tests 

Once the numerical tests have been cleaned up, the right-hand side should always be computed with (optional) augmented integration. 
There should be a parameter 'r_plus' to control the added interpolation quality of the right-hand side. 
Notably, if 'r_plus == 0', then there should be a fallback that avoid repeated computation of the mass matrix.
Similarly, the errors should be computed with augmented integration.

# (INACTIVE) Speed computation of conjugation in sparse matrices

whitney2D/poissonmixedbc2Da.cpp takes too long to assemble the matrices.
Try out a subroutine to reduce the computational effort. 

# (INACTIVE) Rename basic to 'base' or 'general' or 'common'

Basic has the wrong connotation, it makes more sense to call it 'base' or 'general'.

Make a survey of a few important projects to get a sense of what name you should use for this one. 
That will give you a sense of what you should do.

Examples: base, common, core, general, std
MFEM:     general 
Feelpp:   core
Lifev:    core 
ngsolve:  std 
Fenics:   common
Hermes:   common 
concepts: ...
https://en.wikipedia.org/wiki/List_of_finite_element_software_packages

# (INACTIVE/ARTICLE) Static analyzer 

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

? -fanalyzer-transitivity
-fanalyzer-verbosity=level # default 

TODO: write an email to the mailing list about the static analyzer. 






# TODO

zeige die ersten zeilen aller header dateien an
for datei in ./*/*.hpp; do head -n 2 $datei; done
include guards umbenennen 
  
class for one-dimensional meshes: graphs with no abandoned nodes
  
hash-tabelle 
  
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
  - Feature liste welche erwuenscht ist: near/far future 
  
unit test layout:
  logstream that is reference to std::cout 
  TEST_DECL_MODULE( str );
  TEST_DECL_CLASS ( str );
  TEST_DECL_BASIC (     );
  TEST_DECL_TOPIC ( str );
  TEST_ANNOUNCE();
  
  









# DONE!

# (DONE) clean up DenseMatrix subsystem 

The following modules look reasonable
- [x] simple scalar functions into the class
- [x] readwrite is never used: retire 
- [x] factorizations (Gauss-Jordan, QR, Cholesky -> inverse )
- [x] complicated operations (transpose,determinant,tensorproduct)
- [x] simple solvers 

# (DONE) Remove dead code 

Make dead code alive again or remove it. Search for the following pieces:

grep 'if(false' ./*/*pp
grep 'if( false' ./*/*pp

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Github badge C++ >= 17

https://img.shields.io/badge/C++-00599C.svg?style=for-the-badge&logo=C++&logoColor=white
https://img.shields.io/badge/-c++-black?logo=c%2B%2B&style=social
https://img.shields.io/badge/C++-14-blue
https://img.shields.io/badge/-C++14-yellow?logo=c%2B%2B
https://img.shields.io/badge/-C++14-cyan?logo=c%2B%2B&style=flat-square
https://img.shields.io/badge/-C++14-skyblue?logo=c%2B%2B&style=flat-square
https://img.shields.io/badge/-C++14-deepskyblue?logo=c%2B%2B&style=flat-square

# (DONE) Feature: Unphysical operations 

For testing purposes. Given a broken differential form
(a) put it into canonical form again.
(b) Alternative, bring it into an equivalent form randomly 

# (DONE) Basic:

[x] Survey of float 
[x] Survey/benchmark of factorial/binomial and unit tests 
[x] program with leak, hello world 
[x] benchmark for memory allocation 

# (DONE) Text output 

- Get text operations to solvers and mesh class 
- combinatorics, vector & operator classes: emphasize text over print 
- text: combinatorics, redirect print 
- text: operators, dense, sparse ; redirect print 
- text: mesh, redirect print 
- redirect through entire code 

# (DONE) add complete constructor interfaces 

Apply the rule of six and declare all constructors explicitly even if merely setting them to default. 

# (DONE) Different elementary solvers 
  
Implement solution algorithms for special matrix types:
- diagonal solve 
- left/right triagonal solve 
- unit left/right/ triagonal solve 
- averages between left and right solves 
  
# (DONE) Precisions for solvers and magic numbers 

The linear algebra solvers may work with any type of precision, such as float or double. 
Replace the 'magic numbers' in the library by multiples of the machine epsilon. 
Generally speaking, try to find a good stopping criterion.

grep -E '([0-9]+([eE][-+]?[0-9]+))' ./**/*cpp ./*/*/*pp
grep -E '([-+]?\.[0-9]+([eE][-+]?[0-9]+)?)' ./*/*pp ./*/*/*pp

# (DONE) Dependencies

- [x] All targets (tests and modules) depend on the header files as well. 
- [x] Moreover, every test depends also on the static/dynamic libraries. 

# (DONE) How to abort

Termination is done either via `abort()` or via `throw(0)`, depending on whether exceptions are disabled or not. 

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Separate build and tests targets 

Introduce two different targets, one for building 
the object files and libraries, and the other for 
building the test files 

all: build tests

# (DONE) VTK OUTPUT 

The general philosophy of the VTK module is to treat VTK as an output format alone. 
  
Simple legacy format:
http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
Reading only ASCII

# (DONE) Implement Hodge star operation 

# (DONE) Rewrite composed operators 

Reduce the code complexity of the composed operators using pointers, flags,
and templated constructors. This reduces code complexity. 
Then basically retire the the clone and heir methods.

# (DONE) enable complex coefficients 

Enable complex coefficients in a minimalist fashion:
compose complex operators from known (real) operators as composed operator
This should be implemented as a composed operator.

# (DONE) Retire the communtativity tests 

The old test files contain commutativity tests. Those should be retired.

# (DONE) Retire unnecessary Whitney tests 

Poisson solvers do not need to be tested for Whitney forms. 
It suffices to keep a complicated example with miex boundary conditions.

# (DONE) Clean out legacy alternative tests in the FEM solver files

Introduce a unit test in solverfem for the Darcy-system; the Maxwell is already there.
First, clean out the non-block systems, then the block systems.
Make sure that everything that is deleted has an analogue in the list of solvers.
The following is recommend:
- cpp Mass: CGM
- csr mass: CGM SSOR
- cpp stiff: 
- csr stiff: minres csr or CGM SSOR
- systems: Herzog-Soodhalter mit operator preconditioning

# (DONE) Introduce a custom check script

Introduce a check script which reports common 'errors' in your cpp file,
that is, stuff you consider important for the design of your code.
For example,

- [ ] Replace assert(false) by a project-specific macro
- [ ] Magic floating point constant in code 

```
grep --line-number --recursively --color 'assert(' ./*pp
grep --line-number --recursively --color 'cout' ./*pp
grep --line-number --recursively --color '.*[0-9]' ./*pp
grep --line-number --recursively --color '[0-9]e' ./*pp
```

# (DONE) OpenMP pragmas conditional compilation

Every occurence of 'pragma omp' should be included with a conditional compilation.
This ensures that no compiler warnings about 'unknown pragmas' are issued when you
compile the code with openMP disabled.

# (DONE) Copy assignment operator for mesh classes

Since the copy constructor is defined, there should also be an assignment operator
for the different mesh classes, or it should be deleted explicitly.

# (DONE) Clean out 'cout' references throughout code 

grep 'cout' ./*/*.?pp 
Conduct a clean out of all direct references 
to the standard preinstalled streams throughout
the entire project so that the above call
should only return very few instances in the main code
and otherwise only stuff in test files.

Moreover, consider replacing all the other stuff
by references to clog instead of cout.

# (DONE) Abgleichen der Gitterweiten bei solverfem 

# (DONE) Unit test for condition numbers of single element matrices 

For dimensions 1, 2, and 3
All form degrees and polynomial degrees up to 6
Construct the single element meshes, build mass and stiffness matrices 
compute their inverses (and check their products)
perform QR algorithm to find the eigenvalues

Requires: regular triangle and tetrahedra

# (DONE) Conditional compilation when openMP

Furthermore, if openMP is enabled, then you should compile with an inclusion of thread-safe random number generation. 

Generally speaking, you should replace explicit instances of 'rand' and 'srand' 
by wrapper functions. This makes it easier to switch to different implementations 
throughout whenever that becomes necessary.

For example:
- random_integer();
- seed_random_integer();

# (DONE) Output of solver component 

The solver component prints should all contain the iteration number if possible.
Each print should start with the iterartion number, followed by the message class, and then all other info
RESTARTED
BREAKDOWN
---------
(NOTICE)
WARNING
INTERIM

# (DONE) Argument names in all header files 
    
The function/method declarations in the header files should provide argument names. 
The names should coincide with the ones in the code but that's not necessary. 

Rationale: this improves readability.

# (DONE) Change the include orders 

Go from the most general down to the most specific.
This ensures any overwriting of macros stays local.
Within each grouping, sort alphabetically.

# (DONE) Define and adopt a custom assert macro

There is a function that performs the assert, 
and a macro that delivers the line number and file name
to a function invocation. No further frills.
Use the custom assert macro throughout the project.

# (DONE) Phantom coordinate in 2D mesh output

Add a phantom coordinate coordinate to the output of 2D meshes to plot functions

# (DONE) Unit test descriptions

Update the unit test **descriptions** in every module. They seem to be off in many cases.

# (DONE) Revise logging output 

The logging procedure needs to be reworked.

In particular, switch to an encapsulated approach: all classes should have the ability
to produce logs of themselves. That way, you can isolate the problem 
in just a few methods throughout the code.

Basically, implement the following methods:

- text:        produces a string presentation (no nl)
- print:       outputs the text() into a given stream (with nl)
- << operator: outputs the text() into a given stream (no nl)
- lg:          outputs the text into the log (with nl).
               This function may take a preamble argument

Revert the current design of logging output: there shouldn't be
any automatic newlines. Instead, re-introduce the newlines in the tests
and deactive the automatic newline in the logging object.

# (DONE) Introduce a LOG switch 

Make the logging framework optional by introducing a macro switch 
that enables/disables the logging framework

Then introduce the logging framework throughout the entire code	uniformly.

This requires that the logging interface should be used in the same way
as the entire script for the logging stuff.

# (DONE) guarded element access 

All objects that feature element access via brackets,
either blocky brackets or round brackets,
also feature an additional at-method with the same effective behavior. 
The difference is that the at-methods 
always perform bound checks,
which may not the case for the bracket access methods.

- Enforce the effective behavior
- Enforce the bound check policy.

# (DONE) Remove dead code 

grep 'if(false' ./*/*pp
grep 'if( false' ./*/*pp

# (DONE) 'threshold' should be renamed 'tolerance'

# (DONE) Improve the iterator interface of IndexRange to allow the full scope

# (DONE) Git ID extraction as macro

# (DONE) Improve the file names used for output

While using the program name for the file output is nice, it is better if you can also add an additional prefix.
That way you can separate the output of different subtasks more easily.

# (DONE) Dynamic library dependencies 

The test programs depend on the object files only if static linking is enabled.
Otherwise, they don't.

# (DONE) Output program info 

While using the program name for the file output is nice, it is better if you can also add an additional prefix.
That way you can separate the output of different subtasks more easily.

# (DONE) Combinatorics generate multiindices tests must actually test something

# (DONE) Dense Matrix rewrite, part 1

- [x] simple solvers remain in one file
- [x] reconcile the norms of dense matrices with the norms of vectors 
- [x] matrix tensor product should be part of the functions
- [x] Cholesky, QR, and Gauss-Jordan should be one file 


# (DONE) Reorganize the legacy directory with subdirectories

There already are a few subdirectories. Extend that a little bit.

# (DONE) ND meshes

Move this into playgrounds together with the unit tests.

# (DONE) static library dependency 
# (DONE) Compilation mode with object files

Make that you can compile everything if the executables only take the object file themselves. 
This is a third compilation mode in addition to static and dynamic libraries.
Make sure it all compiles. 

# (DONE) Prevent warnings from the external stb libraries

# (DONE) Go over the manuals of GCC and Clang, add more possible warnings 

Introduce a larger amount of warnings. Only use those that are not enabled by default. 

However Turn off the following warnings:
- [x] -Wc++98-compat-local-type-template-args
- [x] -Wreserved-identifier
- [x] -Wold-style-cast
- [x] -Wcovered-switch-default

Following this, go over the list of warnings and re-order everything for the sake of consistency. 
Check what needs to be retired.

# (DONE) Finish the printing in the nullspace computation

Enable for all nullspace vectors printing for any polynomial degree.

# (DONE) Printing of higher order & Clean up of unit tests 

Most routines only print if r == 1. Generalize that.
- [x] enable higher-order printing wherever convenient, and provide higher-order printing. Agree on polydegree
- [x] writeCellVector data: print barycentric 2-forms 
- [x] writeCellScalar data: print barycentric n-forms
- [x] Apply uniform format to Darcy, Maxwell, and curlcurl
- [x] Apply uniform format to lshaped?
- [x] Apply uniform format to Poisson

# (DONE) speed up canonicalize?

- [x] Reduce the number of zero entries in the canonicalization
- [x] Adapt the randomization accordingly

# (DONE) Fix finite difference tests

Check for notanumber and go over the different tests to ensure they test something.

# (DONE) Robustness under shuffling combinatorial data

The code assumes at several points that indexmaps and sigmas are ordered, 
at least when the indexmaps contains only one element. 

Result: 
 - Shuffling the multindices is fine. 
 - Shuffling the sigmas is more difficult.

# (DONE) shake the coordinates in tests where there is no explicit functions living on them 
- [x] fem 
- [x] solvers?
- [x] several finite element tests

# (DONE) Either implement LQ factorization or retire it completely

Implement the LQ factorization and test it

# (DONE) Some warnings to process:

Understand why a compiler might warn about weak vtables and how to avoid that issue. 
This concerns IndexMap and MultiIndex in particular. 
https://stackoverflow.com/questions/23746941/what-is-the-meaning-of-clangs-wweak-vtables

This is not particulary dangerous for this project.

# (DONE) Averaging for Sullivan and Whitney spaces

For each element, you extract the coefficients of the local basis associated with a subsimplex using a Gram-Matrix. 
This works for Sullivan and Whitney bases alike, no difference. 
You can then average according to some scheme, such as:
- [x] arbitrary choice
- [x] weight uniformly
- [x] weight by volume 
Thus you can always average from the non-conforming into the conforming space.

# (DONE) dependencies for object file compilation  

Ensure that object files have all their dependencies correctly defined.

# (DONE) fix warnings about printf truncation 

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


# (DONE): Choose the variant of the determinant. Make manual computation a special case 

# (DONE): The finite element matrices should allow the case r==0 for volume forms 

# (DONE): unit test for LQ factorization 

# (DONE): Fix output of generate multiindices 


