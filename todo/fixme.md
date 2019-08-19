


# TODO short term and minor 

## Copy assignment operator for mesh classes

Since the copy constructor is defined,
there should also be an assignment operator
for the different mesh classes,
or it should be deleted explicitly.

## Logging class 

Even though advanced logging control would be desirable, 
for the time being it is ufficient if the logging capabilities 
are merely present.

- First layer: semantic wrappers for the cpp streams 
- Second layer: advanced logging classes for the alias streams
- Third layer: primitive MACROS that wrap


## Matrix Market subproject 

This subproject should be extended systematically.


## interesting maps

Use the US states map from Randy's source code 
and implement it here. Try to find other triangulations 
too and integrate them as examples. 


## Debug longest edge bisection 

**DONE**
The current implementation uses a stack.
Due to elements appearing multiple times
and being bisected sequentially multiple times,
it is necessary to replace the stack container 
with a list container (or whatever) and 
simulate the stack behavior manually. 


## rename todo file in matrix market subdirectory 

**DONE**


## General check routine for mesh class 

**DONE**
- [x] bisection of single edges distributed throughout the mesh
- [x] longest vertex bisection 2D
- [x] longest vertex bisection 3D
- [x] NVB 2D
- [x] NVB 3D


## What else?
