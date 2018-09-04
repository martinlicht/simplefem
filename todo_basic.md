  
  
  # General TODO-List for the entire project 
  
  ## define a custom attest macro
    
     There is a function that performs the assert, 
     and a macro that delivers the line number and file name
     to a function invocation. No further frills.
    
  ## Adoption of custom assert macro
    
     Use the custom assert macro throughout the library.
  
  
  ## guarded element access 
  
    All objects that feature element access via brackets,
    either blocky brackets or round brackets,
    also feature an additional at-method with the same effective behavior. 
    The difference is that the at-methods 
    always perform bound checks,
    which may not the case for the bracket access methods.
    
    x Enforce the effective behavior
    x Enforce the bound check policy.
  
  
  ## implement minimalist file stream wrapper 
  
    openinputfile( std::string );
    openoutputfile( std::string );
    
  
  
  ## Argument names in all header files 
    
    The function/method declarations in the header files 
    should provide argument names. 
    The names shall coincide in the header and code files. 
  
  
  
  
  
  
  
  
