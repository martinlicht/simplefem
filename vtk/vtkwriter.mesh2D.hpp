#ifndef INCLUDEGUARD_VTK_MESH2D_WRITER
#define INCLUDEGUARD_VTK_MESH2D_WRITER

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

class VTK_MeshWriter_Mesh2D
{
    
    public:
        
        VTK_MeshWriter_Mesh2D( Mesh& m2D, std::ostream& os, const std::string& name );
        
        VTK_MeshWriter_Mesh2D writeCoordinateBlock();
        VTK_MeshWriter_Mesh2D writeCoordinateBlock( const FloatVector& );
        
        VTK_MeshWriter_Mesh2D writeTopDimensionalCells();
        
        VTK_MeshWriter_Mesh2D writeVertexScalarData(
          const FloatVector&, 
          const char* name, Float scaling = 1.
        );
        
        VTK_MeshWriter_Mesh2D writeCellScalarData(
            const FloatVector&, 
            const char* name, Float scaling = 1. 
        );
        
        VTK_MeshWriter_Mesh2D writeCellVectorData(
            const FloatVector& x, 
            const FloatVector& y, 
            const FloatVector& z, 
            const char* name, Float scaling = 1. 
        );
        
    private:
        
        VTK_MeshWriter_Mesh2D writePreamble( const std::string& name );
    
        const Mesh& mesh;
        std::ostream& os;
        
};

#endif
