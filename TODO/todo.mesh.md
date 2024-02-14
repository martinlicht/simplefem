
# (LOW) Bisection method return value 

The bisection methods should return the number of triangles that have been bisected. 
That should not depend on the order of the edges being listed for marking,
that can essentially just be a set.
Generally speaking, it should be a transversable container. 



# Mesh module

The mesh module provides an abstract base class for the management of meshes
and several implementations 


  - Mesh
    - manifold 2D
    - manifold 3D
    - allgemeine manifold
    - allgemeines mesh
  - uniform refinement:
      zuerst verstehen, dann implementieren
  - bisection refinement:
      zuerst verstehen, dann implementieren
  - Ein eigenes Ausgabeformat ist notwendig
      Unterscheidung nach Spezialklassen
      I/O in andere Formate ist zweitrangig
  
  
  


  o Coordinates 
    - check 
    - robusten 
    
  o Abstract Mesh Class
    - simplify the interface and the access to super simplices 
    - move the coordinates member somewhere different 
    
  o Mesh 2D Triangulation 
    - clean up the naming of variables
      (temporarily copy the class declaration into the cpp file, edit, copy back)
    - develop vertex-parent lists as local data.
    - improve upon the bisection methods
    - introduce a for-each interface to improve the navigation  
  
  
  
# rename Mesh class to MeshInterface or something similar
    
    The general mesh class will serve as a general interface 
    for all finite element methods. The details of implementation 
    are deferred to the particular classes of meshes 
    that provide the data structures, IO, and refinement.
    The common mesh interface will provide:
    - check feature availability
    - navigate according to features 
    - but NOT ANY mesh manipulation
  
  
  o One-dimensional mesh
  
    A class of meshes equivalent to a graph without abandoned nodes. 
    Local refinement and (consecutively) global refinement.
    - implement local and global refinement
  
  
  o mesh input/output:
    
      - coordinates io -> readwrite.coordinates
      - readwrite      -> mesh.manifold2D
      - use your own file format for the time being
    - argument namen auch in den header files
    - argument namen abgleichen 
  
  
  o implement IO for each concrete mesh class 
  
    Each mesh class that contains actual data shall provide 
    an input/output functionality. The data of the mesh 
    are written AS IS into a stream and loaded AS IS
    from a stream. Flags shall control the level of detail 
    that is being written. 
    DEPENDENCY: file stream wrappers 
  
  
  o split off coordinate object 
  
    Coordinates become a fully independent object.
    They are temporarily owned by a mesh which forwards them 
    to applications and in turn conducts manipulations to the 
    vertices in accordance to the changes of its own topology.
  
  
  o Three-dimensional meshes 
  
    A class for three-dimensional manifold-like meshes.
    Implement IO interface.
    Implement global uniform refinement without any optimization
    Implement bisection local method.
  
  
  o uniform refinement 2D
    implement the two dimensional refinement algorithm 
    with optimal complexity.
  
  
  o bisection refinement 2D
    check whether the complexity is optimal in each case.
  
  o longest edge refinement in 2D
  
  o bisection algorithm in 3D
  
  o longest edge refinement in 3D
  
  
  
  o implement submesh class 
    - These objects inherit the mesh interface and implement it. 
    - Their internal function consists of the capabilities.  
      of pointing to another mesh. 
    - pointers to the original mesh objects so that data 
      can be transfered.

  o implement patch constructor 
    - The star constructor will construct submeshes from a given mesh


