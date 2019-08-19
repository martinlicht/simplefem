

#  VTK OUTPUT 

The general philosophy of the VTK module 
is to treat VTK as an output format alone. 


  o improved VTK output
    - output the header information
    - output meshes 
    - output mesh data (vertex values, etc.)
    
  
  
  ************* VTK Module
  - Implement a reader + writer class
    and classes for the different parts of a VTK-File
  - Read/Write VTK Headers
  - Read/Write data sets (derived from generic data set class)
  - Perhaps conversion between data set classes
  - Conversion between non-VTK related formats
  - Direct IO with mesh and algebra classes via specialized routines 
  
  Simple legacy format:
    http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    Reading only ASCII
  
  
