	
	 
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
      
      