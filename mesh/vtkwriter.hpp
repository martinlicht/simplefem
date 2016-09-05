
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
        
        // writeVertexData();
        
    private:
    
        const SimplicialMesh& mesh;
        std::ostream& os;
};