
# MODULE Operators

This module provides interfaces for abstract linear algebra.
The two core classes are LinearOperator and FloatVector.

FloatVector captures the idea of a vector in linear algebra
and is just a fixed-length array of floating-point numbers.

LinearOperator is an abstract base class for linear operators 
that act on vectors. 

The module provides some further instantiations for the linear operator class.


## TODO: For-each interface

Equip FloatVector with a for-each interface.
  [*] Code written and compiles 
  [ ] Write unit tests 


## TODO: abstract float vector interface

*Do in conjunction with the DenseMatrix class*

Introduce FloatVector as a base class without internal memory 
but instead with two lambdas that simulate the element read/write methods.
In fact, you can just leave both as optional arguments. 
For example, so it may be possible to leave the write method undefined.

Finally, introduce a new derived class, say, PhysicalFloatVector,
that uses an actual internal container for the storage. 



## TODO: 





	
	 
    FloatVector
    * initializer list constructor
    * check
    * check im constructor 
    * unittest komplett 
    
    LinearOperator
    * check
    * check in all methods 
    * (no unittest because purely virtual)
    - function wrapper 
    
    ScalingOperator
    * check
    * check im constructor 
    - unittest komplett 
    
    DiagonalOperator
    * check
    * check im constructor 
    - unittest komplett 
    
    SumOperator
    * check
    * check im constructor 
    * unittest with scaling operator
    
    ProductOperator
    * check
    * check im constructor 
    * unittest with scaling operator
    
    -----------------------------------------
    
    FloatVector
    - interne daten umstellen auf C-style arrays 
      C++ vector nur noch als kopie zurückgeben
    - high accuracy scalar product (assembler)
    - high speed scalar product (assembler)
    - getInternalPointer 
    
    DiagonalOperator
    - interne daten umstellen auf C-style arrays 
    
    Sum and Product operators 
    - Einige derivierte operatorklassen haben spezialisierte summen und produkte.
      Um diese kompatibel zu halten mit den allgemeinen Summen und Produkten,
      wollen wir einen Cast-Mechanismus verwenden. Eine Option ist es, 
      intern zu speichern dass eine Summe / ein Produkt als Wrapper fungiert 
      eines anderen Operators 
    
    Allgemein:
    - Move-Semantik und weitere fortgeschrittene Konzepte sollten helfen, 
      den Code zu verbessern und automatische Optimierung greifbarer werden lassen.
      
      
