

/**/

#include <cmath>

#include <fstream>
#include <functional>



#include "../../base/include.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/global.inclsullivan.hpp"
#include "../../fem/utilities.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 2D Maxwell System on L-shaped domain" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = LShapedDomain2D();
        
        M.getCoordinates().scale( 1.0 );
        
        M.check();
        
        M.automatic_dirichlet_flags();

        
        
        std::function<FloatVector(const FloatVector&)> function_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 2 );
                return FloatVector({
                    -1.0
                    ,
                    +1.0
                    });
            };

            
            
            
        

        
        

        const int min_l = 0; 
        
        const int max_l = 4;
        
        const int min_r = 1; 
        
        const int max_r = 1;
        

        
        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ )
        {
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ )
            {
                
                LOG << "Level: "             << min_l << " <= " << l << " <= " << max_l << nl;
                LOG << "Polynomial degree: " << min_r << " <= " << r << " <= " << max_r << nl;
                        
                LOG << "... assemble matrices" << nl;
        
                
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
                SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                
                
                auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                mat_A.sortandcompressentries();
                
                auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
                mat_B.sortandcompressentries();
                
                auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                mat_C.sortandcompressentries();
                
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                auto C  = MatrixCSR( mat_C  );

                
                {

                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r,   function_rhs  );
                    
                    FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );
                    
                    FloatVector sol( vector_incmatrix.getdimin(), 0. );
                    
                    
                    timestamp start = timestampnow();
                    
                        LOG << "... iterative solver" << nl;
                        
                        auto PA = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )
                                    + MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix );
                        auto PC = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                                    + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix );
                        
                        LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                        LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;

                        const auto PAinv = inv(PA,desired_precision,-1);
                        const auto PCinv = inv(PC,desired_precision,-1);

                        FloatVector  x_A( A.getdimin(),  0. ); 
                        FloatVector& x_C = sol;
                        
                        const FloatVector  b_A( A.getdimin(),  0. ); 
                        const FloatVector& b_C = rhs; 
                        
                        BlockHerzogSoodhalterMethod( 
                            x_A, 
                            x_C, 
                            b_A, 
                            b_C, 
                            -A, Bt, B, C, 
                            desired_precision,
                            1,
                            PAinv, PCinv
                        );

                        if(false){ // TODO(martinlicht): fix 
                    
                            FloatVector res = sol;
                            
                            HodgeConjugateResidualSolverCSR_SSOR( 
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),   A.getC(),  A.getV(), 
                                B.getA(),   B.getC(),  B.getV(), 
                                Bt.getA(), Bt.getC(), Bt.getV(), 
                                C.getA(),   C.getC(),  C.getV(), 
                                res.raw(),
                                desired_precision,
                                1,
                                desired_precision,
                                0
                            );
                        
                        }
                    
                    timestamp end = timestampnow();
    
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    assert( sol.is_finite() );

                    LOG << "... compute residual" << nl;

                    Float residualnorm      = ( rhs - B * inv(A,desired_precision) * Bt * sol - C * sol ).norm();

                    LOG << "residual:   " << residualnorm  << nl;

                    
                    auto computed_sol  = vector_incmatrix * sol;
                    
                    {
                        std::fstream fs( get_available_filename(get_basename(__FILE__)), std::fstream::out );
                        VTKWriter vtk( M, fs, get_basename(__FILE__) );

                        vtk.write_cell_vector_data( function_rhs, "rhs_field" );
                        
                        {
                            auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );
                            auto lowest_sol = interpol_matrix * computed_sol;
                            auto lowest_rhs = interpol_matrix * interpol_rhs;
                            vtk.write_cell_vector_data_barycentricgradients( lowest_rhs,   "righthandside" );
                            vtk.write_cell_vector_data_barycentricgradients( lowest_sol,   "solution_calculation" );

                            auto constant_field = [](const FloatVector& vec) -> FloatVector {
                                assert( vec.getdimension() == 2 );
                                return FloatVector({ 1., 1. });
                            };
                            
                            auto interpolated_constant_r = Interpolation( M, M.getinnerdimension(), 1, r, constant_field );
                            auto interpolated_constant_0 = Interpolation( M, M.getinnerdimension(), 1, 0, constant_field );

                            vtk.write_cell_vector_data_barycentricgradients( interpol_matrix * interpolated_constant_r, "constant_r" );
                            vtk.write_cell_vector_data_barycentricgradients(                   interpolated_constant_0, "constant_0" );
                        }

                        {
                            auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r+1 );
                            auto lowest_ndiv = interpol_matrix * scalar_incmatrix * x_A;
                            vtk.write_cell_scalar_data( lowest_ndiv, "negative_divergence" );
                        }
                        
                        {
                            auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r-1 );
                            auto lowest_curl = interpol_matrix * vector_diffmatrix * computed_sol; 
                            vtk.write_cell_scalar_data_barycentricvolumes( lowest_curl, "curl" );
                        }

                        fs.close();
                    }


                }
                
            }

            if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
    
        } 
                
    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}




