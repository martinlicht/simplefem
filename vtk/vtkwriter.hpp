#ifndef INCLUDEGUARD_VTK_VTKWRITER
#define INCLUDEGUARD_VTK_VTKWRITER

#include <functional>
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
        VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const FloatVector& z );
        VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const std::function<Float(int)>& func_z );
        
        
        VTKWriter write_vertex_scalar_data( const std::function<Float(int)>& datafunction,                  const std::string& name, Float scaling = 1. );
        VTKWriter write_vertex_scalar_data( const std::function<Float(const FloatVector&)>& function,       const std::string& name, Float scaling = 1. );
        VTKWriter write_vertex_scalar_data( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling = 1. );
        VTKWriter write_vertex_scalar_data( const FloatVector& pointvalues,                                 const std::string& name, Float scaling = 1. );
        
        VTKWriter write_cell_scalar_data( const std::function<Float(int)>& datafunction,                  const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_scalar_data( const std::function<Float(const FloatVector&)>& function,       const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_scalar_data( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_scalar_data( const FloatVector& cellvalues,                                  const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_scalar_data_barycentricvolumes(
                                          const FloatVector& volumevalues,                                const std::string& name, Float scaling = 1. );
        
        VTKWriter write_cell_vector_data( const std::function<FloatVector(int)>& datafunction,            const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_vector_data( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_vector_data( const FloatVector& x, 
                                          const FloatVector& y, 
                                          const FloatVector& z,                                           const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_vector_data_barycentricgradients(
                                          const FloatVector& gradvalues,                                  const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_vector_data_barycentriccrosses(
                                          const FloatVector& crossvalues,                                 const std::string& name, Float scaling = 1. );
        VTKWriter write_cell_vector_data_Euclidean(
                                          int outerdim, 
                                          const FloatVector& directions,                                  const std::string& name, Float scaling = 1. );

        VTKWriter write_point_cloud( const DenseMatrix& coords );


    private:
        
        VTKWriter write_preamble( const std::string& name );
    
        VTKWriter write_coordinate_block();
        VTKWriter write_coordinate_block( const FloatVector& z );
        VTKWriter write_coordinate_block( const std::function<Float(int)>& func_z );
        
        VTKWriter write_top_dimensional_cells();
        
        
        const Mesh& mesh;
        std::ostream& os;

        enum class Stage : int8_t {
            nothing      = -1,
            preamble     =  0,
            coordinate   =  1,
            cells        =  2,
            fielddata    =  3,
            appendix     =  4 
        };
        
        Stage current_stage;
        
        
};

#endif
