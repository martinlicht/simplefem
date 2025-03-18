

/**/

#include <cmath>

#include <fstream>
#include <functional>


#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../utility/math.hpp"
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
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D Maxwell System" << nl;
    
    if(true){

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.getCoordinates().scale( 1.1 );
        
        M.check();
        
        M.automatic_dirichlet_flags();

        
        
        std::function<FloatVector(const FloatVector&)> experiment_sol = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 3 );
                // return FloatVector({ 1. });
                return FloatVector({ 
                    bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                    ,
                    bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                    ,
                    bumpfunction(vec[0])*bumpfunction(vec[1])*bumpfunction(vec[2])
                });
            };
        

        std::function<FloatVector(const FloatVector&)> experiment_ndiv = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 3 );
                // return FloatVector({ 1. });
                return FloatVector( { 
                    - bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                    - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2])
                    - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2])
                });
            };
        

        // + a_y dyx + a_z dzx + b_x dxy + b_z dzy + c_x dxz + c_y dyz
        // - a_y dxy - a_z dxz + b_x dxy - b_z dyz + c_x dxz + c_y dyz
        // - a_y dxy + b_x dxy - a_z dxz + c_x dxz - b_z dyz + c_y dyz
        // 

        std::function<FloatVector(const FloatVector&)> experiment_curl = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 3 );
                // return FloatVector({ 1. });
                return FloatVector( { // - partial_y + partial_x
                    - bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2]) // xy
                    + bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                    ,
                    - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2]) // xz
                    + bumpfunction_dev(vec[0]) * bumpfunction(vec[1]) * bumpfunction(vec[2])
                    ,
                    - bumpfunction(vec[0]) * bumpfunction(vec[1]) * bumpfunction_dev(vec[2]) // yz
                    + bumpfunction(vec[0]) * bumpfunction_dev(vec[1]) * bumpfunction(vec[2])
                });
            };            

        std::function<FloatVector(const FloatVector&)> experiment_rhs = 
            [=](const FloatVector& vec) -> FloatVector {
                Assert( vec.getdimension() == 3 );
                
                return FloatVector({
                    -
                    bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                    ,
                    -
                    bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                    ,
                    -
                    bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction_devdev(vec[1]) *        bumpfunction(vec[2])
                    -
                    bumpfunction(vec[0])        * bumpfunction(vec[1]) *        bumpfunction_devdev(vec[2])
                    });
            };

            
            
            
        

        ConvergenceTable contable("Mass error and solver residual");
        
        contable << "sigma_error" << "u_error" << "du_error" << "residual" << "time" << nl;
        
        

        const int min_l = 0; 
        
        const int max_l = 3;
        
        const int min_r = 2; 
        
        const int max_r = 2;
        

        
        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ )
        {
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ )
            {
                
                LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                
                LOG << "... assemble mass matrices" << nl;

                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );

                LOG << "... assemble inclusion matrices" << nl;
                
                SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                SparseMatrix pseudo_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();
                
                LOG << "... assemble algebraic matrices" << nl;

                SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                SparseMatrix vector_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 );
                SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                SparseMatrix pseudo_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
                SparseMatrix pseudo_elevationmatrix_t = pseudo_elevationmatrix.getTranspose();

                
                
                
                LOG << "... compose system matrices" << nl;

                auto A  = Conjugation( MatrixCSR(scalar_massmatrix), MatrixCSR(scalar_incmatrix) );

                auto Bt = MatrixCSR(scalar_incmatrix_t) & MatrixCSR(scalar_diffmatrix_t) & MatrixCSR(vector_elevationmatrix_t) & MatrixCSR(vector_massmatrix) & MatrixCSR(vector_incmatrix); // upper right
                
                auto B  = Bt.getTranspose(); //pseudo_incmatrix_t & pseudo_massmatrix & diffmatrix & vector_incmatrix; // lower left
                
                auto C  = Conjugation( MatrixCSR(pseudo_massmatrix), MatrixCSR(pseudo_elevationmatrix) & ( MatrixCSR(vector_diffmatrix) & MatrixCSR(vector_incmatrix) ) );
                
                auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                
                
                
                LOG << "... compose preconditioner data" << nl;
                
                auto PA = Conjugation( MatrixCSR(scalar_massmatrix), MatrixCSR(scalar_incmatrix) )
                          +
                          Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(vector_elevationmatrix) & ( MatrixCSR(scalar_diffmatrix) & MatrixCSR(scalar_incmatrix) ) );
                
                auto PC = Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(vector_incmatrix) )
                          + 
                          Conjugation( MatrixCSR(pseudo_massmatrix), MatrixCSR(pseudo_elevationmatrix) & ( MatrixCSR(vector_diffmatrix) & MatrixCSR(vector_incmatrix) ) );
                
                LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
            
                {

                    LOG << "... interpolate explicit solution and rhs" << nl;
                    
                    const auto& function_ndiv = experiment_ndiv;
                    const auto& function_sol  = experiment_sol;
                    const auto& function_curl = experiment_curl;
                    const auto& function_rhs  = experiment_rhs;
                    
                    FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r, function_ndiv );
                    FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r, function_sol  );
                    FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r, function_curl );
                    
                    FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r, function_rhs  );
                    
                    FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );

                    FloatVector sol( vector_incmatrix.getdimin(), 0. );
                    
                    FloatVector  x_A( A.getdimin(),  0. );  // x_A is a dummy for the gradient 
                    FloatVector& x_C = sol;                 // x_C is a reference to the solution 
                    
                    const FloatVector  b_A( A.getdimin(),  0. ); 
                    const FloatVector& b_C = rhs; 
                    
                    LOG << "... iterative solver" << nl;
                        
                    timestamp start = timestampnow();

                    if(false){
                        auto Z  = MatrixCSR( B.getdimout(), B.getdimout() ); // zero matrix
                        const auto PAinv = inv(PA,desired_precision,-1);
                        const auto PCinv = inv(PC,desired_precision,-1);
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
                    }

                    { // TODO(martinlicht): fix 
                        HodgeConjugateResidualSolverCSR_SSOR( 
                            B.getdimout(), 
                            A.getdimout(), 
                            x_C.raw(), 
                            rhs.raw(), 
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

                    
                    assert( sol.is_finite() );

                    auto ndiv = inv(A,desired_precision) * Bt * sol;
                    
                    auto curl = pseudo_elevationmatrix * vector_diffmatrix * vector_incmatrix * sol;
                    
                    LOG << "... compute error and residual" << nl;

                    auto errornorm_aux_ndiv = interpol_ndiv - scalar_incmatrix * ndiv;
                    auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix *  sol;
                    auto errornorm_aux_curl = interpol_curl -                    curl;

                    Float errornorm_ndiv_sq = ( errornorm_aux_ndiv * ( scalar_massmatrix * errornorm_aux_ndiv ) );
                    Float errornorm_sol_sq  = ( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                    Float errornorm_curl_sq = ( errornorm_aux_curl * ( pseudo_massmatrix * errornorm_aux_curl ) );
                    Float residualnorm      = ( rhs - B * inv(A,desired_precision) * Bt * sol - C * sol ).norm();

                    LOG << errornorm_ndiv_sq << space << errornorm_sol_sq << space << errornorm_curl_sq << nl;

                    assert( errornorm_ndiv_sq >= 0 );
                    assert( errornorm_sol_sq  >= 0 );
                    assert( errornorm_curl_sq >= 0 );

                    Float errornorm_ndiv = std::sqrt( errornorm_ndiv_sq );
                    Float errornorm_sol  = std::sqrt( errornorm_sol_sq  );
                    Float errornorm_curl = std::sqrt( errornorm_curl_sq );
                    
                    LOG << "div  error sq: " << errornorm_ndiv_sq << nl;
                    LOG << "sol  error sq: " << errornorm_sol_sq << nl;
                    LOG << "curl error sq: " << errornorm_curl_sq << nl;
                    
                    LOG << "div  error: " << errornorm_ndiv << nl;
                    LOG << "sol  error: " << errornorm_sol << nl;
                    LOG << "curl error: " << errornorm_curl << nl;
                    
                    LOG << "residual:   " << residualnorm  << nl;

                    contable << errornorm_ndiv;
                    contable << errornorm_sol;
                    contable << errornorm_curl;
                    contable << residualnorm;
                    contable << Float( end - start );
                    contable << nl;

                    contable.lg();

                    

                    {
                        std::fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                        VTKWriter vtk( M, fs, getbasename(__FILE__) );

                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                            const auto printable_sol = interpol_matrix * scalar_incmatrix * ndiv; 
                            vtk.write_cell_scalar_data( printable_sol, "computed_negative_divergence" , 1.0 );
                        }
                    
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );
                            const auto printable_sol = interpol_matrix * vector_incmatrix * sol; 
                            vtk.write_cell_vector_data_barycentricgradients( printable_sol, "computed_solution" , 1.0 );
                        }
                    
                        {
                            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r-1 );
                            const auto printable_curl = interpol_matrix * vector_diffmatrix * vector_incmatrix * sol; 
                            Assert( printable_curl.getdimension() == 6 * M.count_simplices(M.getinnerdimension()), 
                                            M.count_simplices(M.getinnerdimension()) );
                            vtk.write_cell_vector_data_barycentriccrosses( printable_curl, "computed_curl" , 1.0 );
                        }
                        
                        vtk.write_cell_scalar_data( function_ndiv, "function_ndiv" , 1.0 );
                        vtk.write_cell_vector_data( function_sol,  "function_sol"  , 1.0 );
                        vtk.write_cell_vector_data( function_curl, "function_curl" , 1.0 );

                        vtk.write_cell_vector_data( function_rhs,  "function_rhs" , 1.0 );
                    
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
