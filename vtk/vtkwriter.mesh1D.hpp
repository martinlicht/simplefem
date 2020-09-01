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
        
        VTK_MeshWriter_Mesh1D( Mesh& m1D, std::ostream& os, const std::string& name );
        
        VTK_MeshWriter_Mesh1D writeCoordinateBlock();
        VTK_MeshWriter_Mesh1D writeCoordinateBlock( const FloatVector& );
        
        VTK_MeshWriter_Mesh1D writeTopDimensionalCells();
        
        VTK_MeshWriter_Mesh1D writeVertexScalarData(
            const FloatVector&, 
            const char* name, Float scaling = 1. 
        );
        
        VTK_MeshWriter_Mesh1D writeCellScalarData(
            const FloatVector&, 
            const char* name, Float scaling = 1. 
        );
        
        VTK_MeshWriter_Mesh1D writeCellVectorData(
            const FloatVector& x, 
            const FloatVector& y, 
            const FloatVector& z, 
            const char* name, Float scaling = 1. 
        );
        
    private:
        
        VTK_MeshWriter_Mesh1D writePreamble( const std::string& name );
    
        const Mesh& mesh;
        std::ostream& os;
        
};

#endif
