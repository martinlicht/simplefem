#ifndef INCLUDEGUARD_VTK_VTKWRITER
#define INCLUDEGUARD_VTK_VTKWRITER

#include <ostream>
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

class VTKWriter
{
    
    public:
        
        VTKWriter( const Mesh& m, std::ostream& os, const std::string& name );
        
        VTKWriter writeCoordinateBlock();
        VTKWriter writeCoordinateBlock( const FloatVector& );
        
        VTKWriter writeTopDimensionalCells();
        
        VTKWriter writeVertexScalarData( const std::function<Float(int)>& datafunction,            const std::string name, Float scaling = 1. );
        VTKWriter writeVertexScalarData( const std::function<Float(const FloatVector&)>& function, const std::string name, Float scaling = 1. );
        VTKWriter writeVertexScalarData( const FloatVector& pointvalues,                           const std::string name, Float scaling = 1. );
        
        
        VTKWriter writeCellScalarData( const std::function<Float(int)>& datafunction,            const std::string name, Float scaling = 1. );
        VTKWriter writeCellScalarData( const std::function<Float(const FloatVector&)>& function, const std::string name, Float scaling = 1. );
        VTKWriter writeCellScalarData( const FloatVector& cellvalues,                            const std::string name, Float scaling = 1. );
        
        VTKWriter writeCellVectorData( const std::function<FloatVector(int)>& datafunction,            const std::string name, Float scaling = 1. );
        VTKWriter writeCellVectorData( const std::function<FloatVector(const FloatVector&)>& function, const std::string name, Float scaling = 1. );
        VTKWriter writeCellVectorData( const FloatVector& x, 
                                       const FloatVector& y, 
                                       const FloatVector& z,                                           const std::string name, Float scaling = 1. );
        VTKWriter writeCellVectorData_Whitney(
                                       const FloatVector& gradvalues,                                  const std::string name, Float scaling = 1. );

        
    private:
        
        VTKWriter writePreamble( const std::string& name );
    
        const Mesh& mesh;
        std::ostream& os;

        enum class Stage {
            nothing      = -1,
            preamble     = 0,
            coordinate   = 1,
            cells        = 2,
            vertexdata   = 3,
            celldata     = 4 
        };
        
        Stage current_stage;
        
        
};

#endif
