

   # Notes and TODO-List for the finite element matrices
   
   This file describes the management of finite element matrices and operators 
   in this project. Generally speaking, at no point do we assume 
   that the operators are actually represented as matrices. Instead
   we merely presume them to be given as abstract operators. 
   
   Another principle throughout this project is that we distingiush
   between elementwise operations and global operations.
   Typically, the global operations are computed using the local operations
   and the only truly global operation is the mapping of coefficients
   from the global space into the local spaces.
   
   Most operations in finite element exterior calculus can be represented
   in terms of canonical spanning sets, and this is precisely what we will 
   do throughout this project. 
   
   We now give a list of the different operations.
   
   
   
   
   ## Element Mass Operator
   
   Gives the mass operator of a simplex 
   in the canonical spanning set 
   of the Pr-family of spaces 
   
   Required input:
   - coordinates of the simplex as input (indicates dimension)
   - form valency
   - polynomial degree
   
   
   ## Element Differential Matrix 
   
   Describes the exterior derivative 
   in the canonical spanning set 
   of the Pr-family of spaces 
   
   Required input:
   - dimension
   - form valency
   - polynomial degree
   
   
   
   ## Element Trace Matrix 
   
   Describes the trace operation 
   from a subsimplex onto a simplex 
   in the canonical spanning set 
   of the Pr-family of spaces 
   
   Required input:
   - dimension of simplex
   - dimension of subsimplex 
   - index of subsimplex 
   - form valency
   - polynomial degree
   
   
   
   ## Global Pr Mapping
   
   Maps the coeffecients associated to the subsimplices 
   for the Pr family of spaces 
   into the direct sum of the discontinuous finite element space
   
   Required input:
   - mesh (includes dimension)
   - form valency
   - polynomial degree 
   
   
   
   ## Global Pr- Mapping
   
   Maps the coeffecients associated to the subsimplices 
   for the Pr- family of spaces 
   into the direct sum of the discontinuous finite element space
   
   Required input:
   - mesh (includes dimension)
   - form valency
   - polynomial degree 
   
   
   
   ## Global Pr Mask
   
   Maps the spaces with boundary conditions 
   into the global space 
   for the Pr family of spaces 
   
   Required input:
   - mesh (includes dimension)
   - form valency
   - polynomial degree 
   
   
   
   ## Global Pr- Mask
   
   Maps the spaces with boundary conditions 
   into the global space 
   for the Pr- family of spaces 
   
   Required input:
   - mesh (includes dimension)
   - form valency
   - polynomial degree 
   
   
   
   
   
   ## Element Extension Operation 
   
   Describes the extension operation
   of the Pr-family of spaces 
   from a subsimplex onto a simplex 
   in the canonical spanning set 
   of the Pr-family of spaces 
   
   Required input:
   - dimension of simplex
   - dimension of subsimplex 
   - index of subsimplex 
   - form valency
   - polynomial degree
   
   
   
   
   
   
   
   
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
  
     
