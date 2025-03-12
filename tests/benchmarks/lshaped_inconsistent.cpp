

/**/

#include <fstream>

#include "../basic.hpp"
#include "../utility/files.hpp"
#include "../utility/convergencetable.hpp"
#include "../operators/composedoperators.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../sparse/matcsr.hpp"
#include "../mesh/mesh.simplicial2D.hpp"
#include "../mesh/examples2D.hpp"
#include "../vtk/vtkwriter.hpp"
#include "../solver/sparsesolver.hpp"
#include "../solver/iterativesolver.hpp"
#include "../solver/inv.hpp"
#include "../solver/systemsparsesolver.hpp"
#include "../solver/systemsolver.hpp"
#include "../fem/global.veewedgehodge.hpp"
#include "../fem/global.diffmatrix.hpp"
#include "../fem/global.elevation.hpp"
#include "../fem/global.interpol.hpp"
#include "../fem/global.massmatrix.hpp"
#include "../fem/global.sullivanincl.hpp"
#include "../fem/utilities.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 2D Maxwell System" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M0 = LShapedDomain2D();
        MeshSimplicial2D M1 = LShapedDomain2D();
        auto& M = M0;

        // M0.automatic_dirichlet_flags();
        // M1.automatic_dirichlet_flags();


        // if(0)
        for( int e = 0; e < M.count_edges(); e++ )
        {
            if( M.count_edge_triangle_parents(e) == 2 ) continue;
            assert( M.count_edge_triangle_parents(e) == 1 );

            int v0 = M.get_edge_vertex( e, 0 );
            int v1 = M.get_edge_vertex( e, 1 );

            auto cs0 = M.getCoordinates().getvectorclone( v0 );
            auto cs1 = M.getCoordinates().getvectorclone( v1 );

            if( cs0[1] == cs1[1] ) {
                M0.set_flag( 0, v0, SimplexFlag::SimplexFlagDirichlet );
                M0.set_flag( 0, v1, SimplexFlag::SimplexFlagDirichlet );
                M0.set_flag( 1, e,  SimplexFlag::SimplexFlagDirichlet );
                LOG << "x: " << e << space << v0 << space << v1 << nl;
            } else if( cs0[0] == cs1[0] ) {
                M1.set_flag( 0, v0, SimplexFlag::SimplexFlagDirichlet );
                M1.set_flag( 0, v1, SimplexFlag::SimplexFlagDirichlet );
                M1.set_flag( 1, e,  SimplexFlag::SimplexFlagDirichlet );
                LOG << "y: " << e << space << v0 << space << v1 << nl;
            }
                
        }

        

        
        std::function<FloatVector(const FloatVector&)> vectorfield_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    -1.0
                    ,
                    +1.0
                    });
            };

        std::function<FloatVector(const FloatVector&)> vectorfield_x = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    1.0
                    ,
                    0.0
                    });
            };

        std::function<FloatVector(const FloatVector&)> vectorfield_y = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    0.0
                    ,
                    1.0
                    });
            };

        std::function<FloatVector(const FloatVector&)> scalarfield_constantone = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    1.0
                    });
            };

            
        

        const int min_l = 0; 
        
        const int max_l = 5;
        
        const int min_r = 1; 
        
        const int max_r = 1;
        

        
        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ )
        {
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ )
            {
                
                LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                
                LOG << "... assemble matrices" << nl;
                
                
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r-1 );
                
                SparseMatrix scalar_massmatrix_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                SparseMatrix scalar_incmatrix0   = FEECSullivanInclusionMatrix( M0, M0.getinnerdimension(), 0, r );
                SparseMatrix scalar_incmatrix1   = FEECSullivanInclusionMatrix( M1, M1.getinnerdimension(), 0, r );
                SparseMatrix scalar_incmatrix0_t = scalar_incmatrix0.getTranspose();
                SparseMatrix scalar_incmatrix1_t = scalar_incmatrix1.getTranspose();

                SparseMatrix scalar_elevationmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, 0, r );
                SparseMatrix scalar_elevationmatrix_t = scalar_elevationmatrix.getTranspose();

                FloatVector dir0 = Interpolation( M, M.getinnerdimension(), 1, 0, vectorfield_x );
                FloatVector dir1 = Interpolation( M, M.getinnerdimension(), 1, 0, vectorfield_y );

                SparseMatrix contractionmatrix0   = FEECBrokenVeeMatrix( M, 2, 1, r-1, 1, 0, dir0 );
                SparseMatrix contractionmatrix1   = FEECBrokenVeeMatrix( M, 2, 1, r-1, 1, 0, dir1 );
                SparseMatrix contractionmatrix0_t = contractionmatrix0.getTranspose();
                SparseMatrix contractionmatrix1_t = contractionmatrix1.getTranspose();
                
                assert( dir0.is_finite() and dir1.is_finite() );
                assert( contractionmatrix0.is_finite() and contractionmatrix1.is_finite() );

                // const int dimension0 = scalar_incmatrix0.getdimin();
                // const int dimension1 = scalar_incmatrix1.getdimin();

                const auto SystemMatrix1 = 
                    Block2x2Operator( 
                        scalar_incmatrix0_t & scalar_diffmatrix_t & contractionmatrix0_t & scalar_massmatrix & contractionmatrix0 & scalar_diffmatrix & scalar_incmatrix0,
                        scalar_incmatrix0_t & scalar_diffmatrix_t & contractionmatrix0_t & scalar_massmatrix & contractionmatrix1 & scalar_diffmatrix & scalar_incmatrix1,
                        scalar_incmatrix1_t & scalar_diffmatrix_t & contractionmatrix1_t & scalar_massmatrix & contractionmatrix0 & scalar_diffmatrix & scalar_incmatrix0,
                        scalar_incmatrix1_t & scalar_diffmatrix_t & contractionmatrix1_t & scalar_massmatrix & contractionmatrix1 & scalar_diffmatrix & scalar_incmatrix1 
                    );
                
                const auto SystemMatrix2 = 
                    Block2x2Operator( 
                            scalar_incmatrix0_t & scalar_diffmatrix_t & contractionmatrix1_t & scalar_massmatrix & contractionmatrix1 & scalar_diffmatrix & scalar_incmatrix0,
                        -( scalar_incmatrix0_t & scalar_diffmatrix_t & contractionmatrix1_t & scalar_massmatrix & contractionmatrix0 & scalar_diffmatrix & scalar_incmatrix1 ),
                        -( scalar_incmatrix1_t & scalar_diffmatrix_t & contractionmatrix0_t & scalar_massmatrix & contractionmatrix1 & scalar_diffmatrix & scalar_incmatrix0 ),
                            scalar_incmatrix1_t & scalar_diffmatrix_t & contractionmatrix0_t & scalar_massmatrix & contractionmatrix0 & scalar_diffmatrix & scalar_incmatrix1 
                    );
                
                const auto SystemMatrix = SystemMatrix1 + SystemMatrix2;
                
                
                
                
                {

                    FloatVector interpol_rhs  = concatFloatVector( -dir0, dir1 );

                    FloatVector interpol_constantone = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield_constantone );
                    
                    FloatVector rhs = concatFloatVector( 
                            scalar_incmatrix0_t * ( scalar_massmatrix_plus * ( -interpol_constantone ) ),
                            scalar_incmatrix1_t * ( scalar_massmatrix_plus * (  interpol_constantone ) )
                        );

                    FloatVector sol( rhs.getdimension(), 0. );


                    
                    timestamp start = timestampnow();
                    
                    auto solver = MinimumResidualMethod( SystemMatrix, desired_precision, SystemMatrix.getdimin() * 1, 1 );
                    solver.solve( sol, rhs );

                    timestamp end = timestampnow();
    
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    assert( sol.is_finite() );

                    
                    
                    
                    if( true ){
                        std::fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );

                        FloatVector z_values( M.count_vertices() );

                        // for( int v = 0; v < M.count_vertices(); v++ ) {
                        //     Float c0 = M.getCoordinates().getdata( v, 0 );
                        //     Float c1 = M.getCoordinates().getdata( v, 1 );
                        //     Float z  = std::sin( c0 * Constants::twopi ) * std::sin( c0 * Constants::twopi );
                        //     z_values[v] = z;
                        // }


                        VTKWriter vtk( M, fs, getbasename(__FILE__) /*, z_values*/ );
                        // vtk.write_coordinate_block();
                        // vtk.write_top_dimensional_cells();

                        auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );

                        auto x_vertexvalues = sol.getslice( 0, scalar_incmatrix0.getdimin()                               );
                        auto y_vertexvalues = sol.getslice(    scalar_incmatrix0.getdimin(), scalar_incmatrix1.getdimin() );

                        auto x_cellvalues = interpol_matrix * scalar_incmatrix0 * x_vertexvalues;
                        auto y_cellvalues = interpol_matrix * scalar_incmatrix1 * y_vertexvalues;

                        if( r == 1 ) vtk.write_vertex_scalar_data( x_vertexvalues, "solution_x_vertex");
                        if( r == 1 ) vtk.write_vertex_scalar_data( y_vertexvalues, "solution_y_vertex");

                        FloatVector interlaced( 3 * x_cellvalues.getdimension() );
                        for( int i = 0; i < x_cellvalues.getdimension(); i++ ) {
                            interlaced[ 3 * i + 0 ] = x_cellvalues[i];
                            interlaced[ 3 * i + 1 ] = y_cellvalues[i];
                            interlaced[ 3 * i + 2 ] = 0.;
                        }


                        vtk.write_cell_vector_data_Euclidean( 3, interlaced, "solution_cell" );
                        vtk.write_cell_vector_data( x_cellvalues, y_cellvalues, FloatVector( M.count_triangles(), 0. ), "solution_cell" );
                        
                        vtk.write_cell_scalar_data( x_cellvalues, "solution_x_cell" );
                        vtk.write_cell_scalar_data( y_cellvalues, "solution_y_cell" );

                        DenseMatrix triangle_midpoints( M.count_triangles(), 2);

                        for( int t = 0; t < M.count_triangles(); t++ ){
                            auto midpoint = M.get_midpoint( 2, t );
                            triangle_midpoints.setrow( t, midpoint );
                        }

                        vtk.write_point_cloud( triangle_midpoints );
                        
                        fs.close();
                    }


                }
                
            }

            if( l != max_l ) { 
                LOG << "Refinement..." << nl; 
                M0.uniformrefinement(); 
                M1.uniformrefinement(); 
            }

    
        } 
                
    }
    
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}




