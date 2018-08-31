  
  
  
    o Pointer or References
      
      Should the use of references be skipped altogether?
      
      
    o Smart Pointers
    
      Should smart pointers be employed throughout the library
      to make it more robust against user malpractice?
    
    
    o namespaces for the project
      
      package everything into a namespace.
      
    
    o extend unit tests 
      
      Make the unit tests for extensive for each class.
      Test every single functionality for its accuracy.
  
  

  **** TODO ***
  
  
  - zeige die ersten zeilen aller header dateien an
    for datei in ./*/*.hpp; do head -n 2 $datei; done
    include guards umbenennen 
  
  Namespaces
    - introduce a namespace for FEECPP
    - subnamespaces below that level 
  
  - class for one-dimensional meshes: graphs with no abandoned nodes
  
  
  
  
  
  
  
  
  - hash-tabelle 
  
  
  Zwischenziele: 
  - zweidimensionales mesh ausgeben und visualisieren 
  - scalar mass matrix aufstellen 
  - funktion interpolieren 
  - interpolante ausgeben 
  --->  
  
  - solve poisson problem with natural boundary conditions 
  - output of solution data, error measurement  
    - interpolation of functions to P\Lambda: punkte aussuchen, polynome auswerten  
  - randbedingungen
  - 
  
  
  
  
  
  
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
    - Feature liste welche erwünscht ist: near/far future 
    
  unit test layout:
    logstream that is reference to std::cout 
    TEST_DECL_MODULE( str );
    TEST_DECL_CLASS ( str );
    TEST_DECL_BASIC (     );
    TEST_DECL_TOPIC ( str );
    TEST_ANNOUNCE();
    
    
    
  
  ************* General layout: 
  
  - Reduce dense matrix module to the core functions of dense matrices
  - Module for Matrix I/O 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  