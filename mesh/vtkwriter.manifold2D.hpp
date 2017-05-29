#ifndef INCLUDEGUARD_VTK_MANIFOLD2D_WRITER
#define INCLUDEGUARD_VTK_MANIFOLD2D_WRITER

#include <iostream>
#include <string>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/manifold.2D.hpp"

/*******************
****  
****  Write a mesh as VTK File 
****  - write Preamble, Coordinate Block, Top-dimensional cells
****  
*******************/

class VTK_MeshWriter_Manifold2D
{
    
    public:
        
        VTK_MeshWriter_Manifold2D( ManifoldTriangulation2D& m2d, std::ostream& os );
        
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
    
        const ManifoldTriangulation2D& mesh;
        std::ostream& os;
        
};

#endif