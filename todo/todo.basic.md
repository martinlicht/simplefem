
# GENERAL TODO LIST FOR THE ENTIRE PROJECT 


## 0. Fixed-size dynamic array and adoption

Define a template class for a dynamically allocated array
whose size cannot be changed after allocation. 
Copy the std::vector interface but do not provide 
resizing and capacity information.

Use that fixed-size array throughout your code whenever appropiate,
replacing the old std::vector variables with the new ones.
This applies in particular to the linear algebra classes.

## 1. define a custom attest macro
    
There is a function that performs the assert, 
and a macro that delivers the line number and file name
to a function invocation. No further frills.



## 2. Adoption of custom assert macro
    
Use the custom assert macro throughout the library.




## 3. guarded element access 

All objects that feature element access via brackets,
either blocky brackets or round brackets,
also feature an additional at-method with the same effective behavior. 
The difference is that the at-methods 
always perform bound checks,
which may not the case for the bracket access methods.

    - Enforce the effective behavior
    - Enforce the bound check policy.


    
##  4. implement minimalist file stream wrapper 

    openinputfile( std::string );
    openoutputfile( std::string );


    
##  5. Argument names in all header files 
    
The function/method declarations in the header files 
should provide argument names. 
The names shall coincide in the header and code files. 

