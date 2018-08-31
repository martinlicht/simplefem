

    General remarks:
    
    NB:
    Verwenden wir spanning sets für die eigentliche Räume
    oder basen für symbolische räume?
    
    We use ordered spanning sets for all vector spaces. 
    That is precisely what defines the matrix structure 
    for any given operator. 

    Element matrices 
    - independent of the dimension
    - use the canonical ordered spanning sets
    - ordering dependent on combinatorial algorithms.

    
    o Construct mass matrix for barycentric polynomials
      in canonical basis.
      
    
    o Construct the barycentric differential matrix 
      The construction requires ordered coordinates of a simplex.
      Should be given rather as the difference matrix
      of the edges, so that it's dimension-independent.
    
    
    o Construct the generalized cofactor matrix
      Necessary for the higher order differential forms.
    
    
    o Construct mass matrix for polynomial differential forms 
      This is only a kronecker product of the former.
      
    
    
    o degree elevation 
      increases the polynomial degree formally 
      
    
    o extension operators 
      - whitney form mit randbedingungen 
      - polynomial sullivan mit randbedingungen
      - erweiterungsoperator
      REMARK:
      As a preparation, understand the extension operators 
      
    
    
    
    Global matrices:
    - rely on a mesh to be constructed 
    - implicitly use an ordered spanners whose ordering depends on
      1. conventions (e.g., order of simplex dimensions) 
      2. the mesh (e.g., order of simplices)
      3. element spanners
    
    
    o Mass matrix 
    o Piecewise derivative 
    o Weighted traces (i.e., jump terms)
    o Average traces (i.e., algebraic traces)
    o For each element matrix the corresponding global matrix
    o Conforming -> non-conforming
    
    
    
    
    
    
    
    ************* FE Matrices: 
      - Static condensation kann einfacher implementiert werden, 
      falls die entsprechenden Freiheitsgrade schon vorher aussortiert werden. 
  
      