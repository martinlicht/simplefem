

/**/

#include <cmath>

#include <fstream>
#include <functional>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.inclsullivan.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D Darcy problem" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();
        M.check_dirichlet_flags();

        
        LOG << "Prepare scalar fields for testing..." << nl;
        

        // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
        
        const Float xfeq = 1.;
        const Float yfeq = 2.;
        const Float zfeq = 3.;
        
        
        // u dx + v dy -> u_y dydx + v_x dxdy = ( v_x - u_y ) dxdy
        
        // phi -> ( - phi_y, phi_x ) -> ( - phi_xx - phi_yy ) dxdy
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] ) });
                // return FloatVector({ blob( vec[0] ) * blob( vec[1] ) * blob( vec[2] ) });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_grad = 
            [=](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                return FloatVector( { 
                        -zfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::cos( zfeq * Constants::twopi * vec[2] ), //xy
                        +yfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::cos( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] ), //xz
                        -xfeq * Constants::twopi * std::cos( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )  //yz
                        // - blob( vec[0] )     * blob( vec[1] )     * blob_dev( vec[2] ),
                        // + blob( vec[0] )     * blob_dev( vec[1] ) * blob( vec[2] ),
                        // - blob_dev( vec[0] ) * blob( vec[1] )     * blob( vec[2] )
                    });
            };
        
        
        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 3 );
                return FloatVector({ 
                    xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    +
                    yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    +
                    zfeq*zfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) * std::sin( zfeq * Constants::twopi * vec[2] )
                    // -
                    //  blob_devdev( vec[0] ) * blob( vec[1] ) * blob( vec[2] )
                    // -
                    // blob( vec[0] ) * blob_devdev( vec[1] ) * blob( vec[2] )
                    // -
                    // blob( vec[0] ) * blob( vec[1] ) * blob_devdev( vec[2] )
                    });
            };
        
        
        

        

        const int min_l = 0; 
        const int max_l = 3;
        
        const int min_r = 1;
        const int max_r = 1;
        
        
        ConvergenceTable contable("Mass error");
        
        contable << "sigma_error" << "u_error" << "sigma_res" << "u_res" << "time" << nl;
        

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
            
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ ) 
            {
                
                LOG << "Polynomial degree: " <<  min_r << " <= " << r << " <= " << max_r << nl;
                    
                LOG << "... assemble mass matrices" << nl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r   );
                
                SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r   );
                
                SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 );

                LOG << "... assemble inclusion matrices" << nl;

                SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r   );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 3, r-1 );
                SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                LOG << "... assemble algebraic matrices" << nl;

                SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                
                LOG << "... compose system matrices" << nl;

                auto A  = Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(vector_incmatrix) );
                
                auto Bt = MatrixCSR(vector_incmatrix_t) & MatrixCSR(diffmatrix_t) & MatrixCSR(volume_massmatrix) & MatrixCSR(volume_incmatrix);
                // upper right
                
                auto B  = Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
                
                auto negA = - A; 
                
                auto C  = MatrixCSR( B.getdimout(), B.getdimout() ); // zero matrix
                
                
                LOG << "... compose preconditioner data" << nl;
                        
                auto PA = A // Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(vector_incmatrix) )
                          +
                          Conjugation( MatrixCSR(volume_massmatrix), MatrixCSR(diffmatrix) & MatrixCSR(vector_incmatrix) );
                
                auto PC = Conjugation( MatrixCSR(volume_massmatrix), MatrixCSR(volume_incmatrix) );

                LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                
                const auto PAinv = inv(PA,desired_precision,-1);
                const auto PCinv = inv(PC,desired_precision,-1);

                {

                    LOG << "... interpolate explicit solution and rhs" << nl;
                    
                    const auto& function_sol  = experiment_sol;
                    const auto& function_grad = experiment_grad;
                    const auto& function_rhs  = experiment_rhs;
                    
                    FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 2, r,   function_grad );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 3, r-1, function_sol  );
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 3, r-1, function_rhs  );
                    
                    FloatVector rhs = volume_incmatrix_t * ( volume_massmatrix * interpol_rhs );

                    FloatVector sol( volume_incmatrix.getdimin(), 0. );
                    
                    FloatVector  x_A( A.getdimin(),  0. );  // x_A is a dummy for the gradient 
                    FloatVector& x_C = sol;                 // x_C is a reference to the solution 
                    
                    const FloatVector  b_A( A.getdimin(),  0. ); 
                    const FloatVector& b_C = rhs; 
                    
                    LOG << "... iterative solver" << nl;
                        
                    timestamp start = timestampnow();

                    {
                        BlockHerzogSoodhalterMethod( 
                            x_A, 
                            x_C, 
                            b_A, 
                            b_C, 
                            negA, Bt, B, C, 
                            desired_precision,
                            1,
                            PAinv, PCinv
                        );
                    }

                    if(false){ // TODO(martinlicht): fix 
                        HodgeConjugateResidualSolverCSR_SSOR( 
                            B.getdimout(), 
                            A.getdimout(), 
                            x_C.raw(), 
                            b_C.raw(), 
                            A.getA(),   A.getC(),  A.getV(), 
                            B.getA(),   B.getC(),  B.getV(), 
                            Bt.getA(), Bt.getC(), Bt.getV(), 
                            C.getA(),   C.getC(),  C.getV(), 
                            FloatVector( C.getdimin(),  0. ).raw(),
                            desired_precision,
                            1,
                            desired_precision,
                            0
                        );
                        x_A = inv(A,desired_precision) * Bt * x_C;
                    }

                    timestamp end = timestampnow();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                                        
                    auto grad = x_A; // inv(A,desired_precision) * Bt * sol;

                    LOG << "... compute error and residual" << nl;

                    auto errornorm_aux_sol  = interpol_sol  - volume_incmatrix *  sol;
                    auto errornorm_aux_grad = interpol_grad - vector_incmatrix * grad;

                    Float errornorm_sol  = std::sqrt( errornorm_aux_sol  * ( volume_massmatrix *  errornorm_aux_sol ) );
                    Float errornorm_grad = std::sqrt( errornorm_aux_grad * ( vector_massmatrix * errornorm_aux_grad ) );
                    Float residual_sol   = ( rhs - B * grad ).norm();
                    Float residual_grad  = ( - A * grad + Bt * sol ).norm();

                    LOG << "error:      " << errornorm_sol << nl;
                    LOG << "grad error: " << errornorm_grad << nl;
                    LOG << "residual:   " << residual_sol << nl;
                    LOG << "residual:   " << residual_grad << nl;

                    contable << errornorm_grad;
                    contable << errornorm_sol;
                    contable << residual_grad;
                    contable << residual_sol;
                    contable << Float( end - start );
                    contable << nl;

                    contable.lg();

                    

                    if( r == min_r ) 
                    {
                        std::fstream fs( get_available_filename(get_basename(__FILE__)), std::fstream::out );
                        VTKWriter vtk( M, fs, get_basename(__FILE__) );

                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r );
                            const auto printable_grad = interpol_matrix * vector_incmatrix * grad; 
                            vtk.write_cell_vector_data_barycentriccrosses( printable_grad, "computed_grad" , 1.0 );
                        }
                    
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 3, 0, r-1 );
                            const auto printable_sol = interpol_matrix * volume_incmatrix * sol; 
                            Assert( printable_sol.getdimension() == (M.getinnerdimension()+1) * M.count_simplices(M.getinnerdimension()), 
                                            printable_sol.getdimension(), M.count_simplices(M.getinnerdimension()) );
                            vtk.write_cell_scalar_data_barycentricvolumes( printable_sol, "computed_sol" , 1.0 );
                        }
                        
                        vtk.write_cell_vector_data( function_grad, "function_grad", 1.0 );
                        vtk.write_cell_scalar_data( function_sol,  "function_sol" , 1.0 );

                        vtk.write_cell_scalar_data( function_rhs,  "function_rhs" , 1.0 );
                    
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
