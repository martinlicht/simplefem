#ifndef INCLUDEGUARD_VTK_MESH3D_WRITER
#define INCLUDEGUARD_VTK_MESH3D_WRITER

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

class VTK_MeshWriter_Mesh3D
{
    
    public:
        
        VTK_MeshWriter_Mesh3D( Mesh& m3D, std::ostream& os, const std::string& name );
        
        VTK_MeshWriter_Mesh3D writeCoordinateBlock();
        
        VTK_MeshWriter_Mesh3D writeTopDimensionalCells();
        
        VTK_MeshWriter_Mesh3D writeVertexScalarData(
          const FloatVector&, 
          const char* name, Float scaling = 1.
        );
        
        VTK_MeshWriter_Mesh3D writeCellScalarData(
            const FloatVector&, 
            const char* name, Float scaling = 1. 
        );
        
        VTK_MeshWriter_Mesh3D writeCellVectorData(
            const FloatVector& x, 
            const FloatVector& y, 
            const FloatVector& z, 
            const char* name, Float scaling = 1. 
        );
        
    private:
        
        VTK_MeshWriter_Mesh3D writePreamble( const std::string& name );
    
        const Mesh& mesh;
        std::ostream& os;
        
};

#endif
