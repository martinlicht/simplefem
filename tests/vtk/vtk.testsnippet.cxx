
#include "../../utility/utility.hpp"
#include "../../operators/floatvector.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../vtk/vtkwriter.hpp"


    
    
inline void internal_print( const Mesh& M, std::string meshname, std::string filename = "" )
{
    
    int inner_dim = M.getinnerdimension();
    
    std::string basename = ( filename == "" ) ? get_basename(__FILE__) : filename;
    std::string usedname = get_available_filename( basename );
    LOG << "Writing into: " << usedname << nl;
    
    std::fstream fs( usedname, std::fstream::out );
    
    VTKWriter vtk( M, fs, meshname );
    // vtk.write_coordinate_block();
    // vtk.write_top_dimensional_cells();
    
    {
        FloatVector V( M.count_simplices(0), 
                        [&M](int i)->Float{
                        FloatVector point = M.getCoordinates().getdata_by_vertex(i);
                        // return i / (Float) M.count_simplices(0);
                        Float x = point[0];
                        Float y = ( ( point.getdimension() >= 2 ) ? point[1] : 1. );
                        Float z = ( ( point.getdimension() >= 3 ) ? point[2] : 1. );
                        return std::sin( Constants::pi * x ) * y + z;
                    });
        
        vtk.write_vertex_scalar_data( V, "vertex_scalar_uno", 1.0 );
    }

    // vtk.write_vertex_scalar_data( 
    //     FloatVector( M.count_simplices(0),  1.5 ),
    //     FloatVector( M.count_simplices(0),  2.5 ),
    //     FloatVector( M.count_simplices(0), -3.5 ),
    //     "vertex_vector_test"
    // );

    vtk.write_cell_scalar_data( FloatVector(M.count_simplices(inner_dim), 2.5 ), "cell_scalar_uno" );
    
    vtk.write_cell_scalar_data( FloatVector(M.count_simplices(inner_dim), [](int){ return random_uniform(); } ), "cell_scalar_duo" );
    
    vtk.write_cell_vector_data( 
        FloatVector( M.count_simplices(inner_dim),  1.5 ),
        FloatVector( M.count_simplices(inner_dim),  2.5 ),
        FloatVector( M.count_simplices(inner_dim), -3.5 ),
        "cell_vector_uno"
    );

    vtk.write_cell_vector_data_Euclidean( 3, FloatVector( M.count_simplices(inner_dim) * 3, [](int i) { 
        switch(i % 3) {
            case 0: return  1.;
            case 1: return -2.;
            case 2: return  3.;
        }
    } ), "cell_vector_duo" );

    fs.close();
}
