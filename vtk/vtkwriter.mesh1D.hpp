#ifndef INCLUDEGUARD_VTK_MESH1D_WRITER
#define INCLUDEGUARD_VTK_MESH1D_WRITER

#include <iostream>
#include <string>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.hpp"

/*******************
****  
****  Write a mesh as VTK File 
****  - write Preamble, Coordinate Block, Top-dimensional cells
****  
*******************/

class VTK_MeshWriter_Mesh1D
{
    
    public:
        
        VTK_MeshWriter_Mesh1D( Mesh& m1D, std::ostream& os );
        
        void writePreamble( const std::string& name );
        void writePreamble( const char* name );
        
        void writeCoordinateBlock();
        
        void writeTopDimensionalCells();
        
        void writeVertexScalarData
        ( const FloatVector&, 
          const char* name, Float scaling = 1. );
        
        void writeCellScalarData
        ( const FloatVector&, 
          const char* name, Float scaling = 1. );
        
        void writeCellVectorData
        ( const FloatVector& x, 
          const FloatVector& y, 
          const FloatVector& z, 
          const char* name, Float scaling = 1. );
        
    private:
    
        const Mesh& mesh;
        std::ostream& os;
        
};

#endif
