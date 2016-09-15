
// #include <cassert>

#include <iostream>
#include <string>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/simplicialmesh.hpp"


/*******************
****  
****  
****  Write a mesh as VTK File 
****  
****  - write Preamble, Coordinate Block, Top-dimensional cells
****  
****  
****  
*******************/



class VTK_MeshWriter
{
    
    public:
        
        VTK_MeshWriter( SimplicialMesh& sm, std::ostream& os );
        
        void writePreamble( const char* name );
        
        void writeCoordinateBlock();
        
        void writeTopDimensionalCells();
        
        void writeVertexScalarData
        ( const FloatVector&, 
          const char* name, Float scaling = 1. );
        
        void writeCellVectorData
        ( const FloatVector& x, 
          const FloatVector& y, 
          const FloatVector& z, 
          const char* name, Float scaling = 1. );
        
    private:
    
        const SimplicialMesh& mesh;
        std::ostream& os;
};