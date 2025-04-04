

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
// #include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Solution of Darcy Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
//             M.automatic_dirichlet_flags();
//             M.check_dirichlet_flags();

            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            // std::function<FloatVector(const std::function<FloatVector(const FloatVector&) ) >scalarfield = 
            
            Float xfeq = 1.;
            Float yfeq = 1.;
            
            
            // u dx + v dy -> u_y dydx + v_x dxdy = ( v_x - u_y ) dxdy
            
            // phi -> ( - phi_y, phi_x ) -> ( - phi_xx - phi_yy ) dxdy
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) });
                };
            
            
            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                             yfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::cos( yfeq * Constants::twopi * vec[1] ), 
                            -xfeq * Constants::twopi * std::cos( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ),
                        });
                };
            
            
            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                        +
                        yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                     });
                };
            
            
            

            

            cout << "Solving Poisson Problem with Neumann boundary conditions" << endl;

            int min_l = 0; 
            int max_l = 7;
            
            int min_r = 3;
            int max_r = 3;
            
            ConvergenceTable contable;
            
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                if( l != 0 )
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    cout << "... assemble matrices" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                    SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();

                    auto mat_A  = vector_incmatrix_t & vector_massmatrix & vector_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = vector_incmatrix_t & diffmatrix_t & volume_massmatrix & volume_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    
                    
                    auto Schur = B * inv(A,1e-10) * Bt;
                    
                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        cout << "...interpolate explicit solution and rhs" << endl;
                        
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r,   function_grad );
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_sol  );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 2, r-1, function_rhs  );
                        
                        cout << "...measure interpolation commutativity" << endl;
            
                        {
                            auto commutatorerror_aux = interpol_rhs - diffmatrix * interpol_grad;
                            Float commutatorerror  = commutatorerror_aux * ( volume_massmatrix * commutatorerror_aux );
                            cout << "algebraic commutator error 1: " << commutatorerror << endl;// << space << commutatorerror2
                        }
                        
//                         {
//                             auto commutatorerror_aux = interpol_grad - vector_massmatrix_inv * diffmatrix_t * volume_massmatrix * interpol_rhs;
//                             Float commutatorerror  = commutatorerror_aux * ( commutatorerror_aux );
//                             cout << "algebraic commutator error 2: " << commutatorerror << endl;// << space << commutatorerror2
//                         }
                        

                        {
                            
                            FloatVector rhs = volume_incmatrix_t * ( volume_massmatrix * interpol_rhs );

                            FloatVector sol( volume_incmatrix.getdimin(), 0. );
                            
                            cout << "...iterative solver" << endl;
                            
                            
                            sol.zero();
                            
                            FloatVector res = sol;
                            


                            timestamp start = gettimestamp();

                            HodgeConjugateResidualSolverCSR_SSOR( // TODO
//                             HodgeConjugateResidualSolverCSR_textbook( 
                                B.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                A.getA(),   A.getC(),  A.getV(), 
                                B.getA(),   B.getC(),  B.getV(), 
                                Bt.getA(), Bt.getC(), Bt.getV(), 
                                res.raw(),
                                1e-10,
                                1
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            
                            
                            auto grad = inv(A,1e-14) * Bt * sol;

                            cout << "...compute error and residual:" << endl;

                            auto errornorm_aux_sol  = interpol_sol  - volume_incmatrix *  sol;
                            auto errornorm_aux_grad = interpol_grad - vector_incmatrix * grad;

                            Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( volume_massmatrix *  errornorm_aux_sol ) );
                            Float errornorm_grad = sqrt( errornorm_aux_grad * ( vector_massmatrix * errornorm_aux_grad ) );
                            Float residualnorm   = ( rhs - B * inv(A,1e-10) * Bt * sol ).norm();

                            cout << "error:     " << errornorm_sol  << endl;
                            cout << "aux error: " << errornorm_grad << endl;
                            cout << "residual:  " << residualnorm  << endl;

                            contable << errornorm_sol;
                            contable << errornorm_grad;
                            contable << nl;

                            contable.print( std::cout );
                            
                        }


                    }
                    
                }

                cout << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();
                
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
