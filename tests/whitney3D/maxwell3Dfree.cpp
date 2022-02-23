

/**/

#include <ostream>
#include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: 3D Maxwell System" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.getcoordinates().scale( 1.1 );
            
            M.check();
            
            M.automatic_dirichlet_flags();

            
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
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
                [=](const FloatVector& vec) -> FloatVector{
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
                [=](const FloatVector& vec) -> FloatVector{
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
                [=](const FloatVector& vec) -> FloatVector{
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
            
            contable << "sigma_error" << "u_error" << "du_error" << "residual";
            
            

            const int min_l = 0; 
            
            const int max_l = 4;
            
            const int min_r = 2; 
            
            const int max_r = 2;
            

            
            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << std::endl;
                LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << std::endl;
                    
                    LOG << "... assemble matrices" << endl;
            
                    
                     // TODO: correct the degrees, perhaps via degree elevation
                    
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );

                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix vector_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 );
                    SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                    SparseMatrix pseudo_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r );
                    SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();
                    
                    SparseMatrix pseudo_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
                    SparseMatrix pseudo_elevationmatrix_t = pseudo_elevationmatrix.getTranspose();

                

                    auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                    auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //pseudo_incmatrix_t & pseudo_massmatrix & vector_elevationmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & pseudo_elevationmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix;
                    mat_C.sortandcompressentries();
                    
                    LOG << "share zero A = " << mat_A.getnumberofzeroentries() << "/" << (Float) mat_A.getnumberofentries() << nl;
                    LOG << "share zero B = " << mat_B.getnumberofzeroentries() << "/" << (Float) mat_B.getnumberofentries() << nl;
                    LOG << "share zero C = " << mat_C.getnumberofzeroentries() << "/" << (Float) mat_C.getnumberofentries() << nl;
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    auto C  = MatrixCSR( mat_C  );

                    auto negA  = A;  negA.scale(-1);
                    auto negB  = B;  negB.scale(-1);
                    auto negBt = Bt; negBt.scale(-1);
                    
                    auto SystemMatrix = C + B * inv(A,desired_precision) * Bt;
                    
                    
                    
                    const auto& foo = inv;
                    
                    
                    
                    {

                        const auto& function_ndiv  = experiment_ndiv;
                        const auto& function_sol  = experiment_sol;
                        const auto& function_curl = experiment_curl;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << endl;
                        
                        FloatVector interpol_ndiv = Interpolation( M, M.getinnerdimension(), 0, r, function_ndiv );
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r, function_sol  );
                        FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r, function_curl );
                        
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r, function_rhs  );
                        
                        FloatVector rhs = vector_incmatrix_t * ( vector_massmatrix * interpol_rhs );

                            
                        FloatVector sol( vector_incmatrix.getdimin(), 0. );
                        
                        
                        
                        
                            


                        if(false)
                        {
                        
                            LOG << "...measure interpolation commutativity" << endl;
                            
//                             auto  commutatorerror_1_aux = interpol_rhs - scalar_diffmatrix * interpol_ndiv - vector_diffmatrix_t * pseudo_massmatrix * interpol_curl;
                            auto  commutatorerror_1_aux
                            = 
                            interpol_rhs
                            - scalar_diffmatrix   * inv(scalar_massmatrix,1e-14) * scalar_diffmatrix_t * vector_massmatrix * interpol_sol
                            - vector_diffmatrix_t * pseudo_elevationmatrix_t * pseudo_massmatrix * pseudo_elevationmatrix * vector_diffmatrix   * interpol_sol;
                            Float commutatorerror_1     = commutatorerror_1_aux * ( vector_massmatrix * commutatorerror_1_aux );
                            LOG << "algebraic commutator error 1: " << commutatorerror_1 << endl;
                            
                            auto  commutatorerror_2_aux = interpol_curl - pseudo_elevationmatrix * vector_diffmatrix * interpol_sol;
                            Float commutatorerror_2     = commutatorerror_2_aux * ( pseudo_massmatrix * commutatorerror_2_aux );
                            LOG << "algebraic commutator error 2: " << commutatorerror_2 << endl;
                            
                            auto  commutatorerror_3_aux = scalar_massmatrix * interpol_ndiv - scalar_diffmatrix_t * vector_elevationmatrix_t * interpol_sol;
                            Float commutatorerror_3     = commutatorerror_3_aux * ( scalar_massmatrix * commutatorerror_3_aux );
                            LOG << "algebraic commutator error 3: " << commutatorerror_3 << endl;
                            
                        }
                        
                        
                        if(false)
                        {
                            
                            LOG << "...iterative solver" << endl;
                            
                            sol.zero();
                            
                            auto X = B * inv(A,1e-14) * Bt + C;

//                             HerzogSoodhalterMethod Solver( C );
                            ConjugateResidualMethod Solver( X );
                            Solver.threshold           = 1e-10;
                            Solver.print_modulo        = 100;
                            Solver.max_iteration_count = 4 * sol.getdimension();

                            timestamp start = gettimestamp();
//                             Solver.solve_fast( sol, rhs );
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << "...compute error and residual:" << endl;

                        }


                        
                        
                        {
                            
                            sol.zero();
                            
//                             rhs = vector_incmatrix_t * vector_diffmatrix_t * pseudo_massmatrix * vector_diffmatrix * interpol_sol;
//                             rhs = B * inv(A,10e-14) * scalar_incmatrix_t * scalar_diffmatrix_t * vector_massmatrix * interpol_sol;
//                             rhs = B * scalar_incmatrix_t * scalar_massmatrix * interpol_ndiv;
                        
                            FloatVector res = rhs;

                            LOG << "...iterative solver" << endl;
                            
                            timestamp start = gettimestamp();

                            LOG << "- mixed system solver" << endl;
//                             if(false)
                            //HodgeConjugateResidualSolverCSR(
                            HodgeConjugateResidualSolverCSR_SSOR(
                            //HodgeConjugateResidualSolverCSR_textbook( 
                                negB.getdimout(), 
                                A.getdimout(), 
                                sol.raw(), 
                                rhs.raw(), 
                                   A.getA(),    A.getC(),    A.getV(), 
                                   B.getA(),    B.getC(),    B.getV(), 
                                  Bt.getA(),   Bt.getC(),   Bt.getV(), 
                                   C.getA(),    C.getC(),    C.getV(),
                                res.raw(),
                                1e-10,
                                100,
                                1e-14, //desired_precision,
                                -1
                            );

                            
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            
                        }

                        if(false)
                        {
                            
                            auto X = Block2x2Operator( A.getdimout() + B.getdimout(), A.getdimin() + Bt.getdimin(), -A, Bt, B, C );

                            FloatVector sol_whole( A.getdimin()  + Bt.getdimin(),  0. );
                            FloatVector rhs_whole( A.getdimout() +  B.getdimout(), 0. );
                            
                            sol_whole.zero();
                            rhs_whole.setslice( 0, A.getdimout(), 0. );
                            rhs_whole.setslice( A.getdimout(), rhs );
                            
                            MinimumResidualMethod Solver( X );
                            Solver.threshold           = 1e-10;
                            Solver.print_modulo        = 500;
                            Solver.max_iteration_count = 10 * sol_whole.getdimension();

                            timestamp start = gettimestamp();
                            Solver.solve( sol_whole, rhs_whole );
                            timestamp end = gettimestamp();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << "...compute error and residual:" << endl;

                            Float residualnorm  = ( rhs_whole - X * sol_whole ).norm();

                            LOG << "combined system residual:  " << residualnorm << endl;

                            sol = sol_whole.getslice( A.getdimout(), B.getdimout() );
                        }


                        assert( sol.isfinite() );

                        auto ndiv = inv(A,1e-14) * Bt * sol;
                        
                        auto curl = pseudo_elevationmatrix * vector_diffmatrix * vector_incmatrix * sol;
                        
                        LOG << "...compute error and residual:" << endl;

                        auto errornorm_aux_ndiv = interpol_ndiv - scalar_incmatrix * ndiv;
                        auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix *  sol;
                        auto errornorm_aux_curl = interpol_curl -                    curl;

                        Float errornorm_ndiv_sq = ( errornorm_aux_ndiv * ( scalar_massmatrix * errornorm_aux_ndiv ) );
                        Float errornorm_sol_sq  = ( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol  ) );
                        Float errornorm_curl_sq = ( errornorm_aux_curl * ( pseudo_massmatrix * errornorm_aux_curl ) );
                        Float residualnorm      = ( rhs - B * inv(A,1e-14) * Bt * sol - C * sol ).norm();

                        LOG << errornorm_ndiv_sq << space << errornorm_sol_sq << space << errornorm_curl_sq << endl;

                        assert( errornorm_ndiv_sq >= 0 );
                        assert( errornorm_sol_sq  >= 0 );
                        assert( errornorm_curl_sq >= 0 );

                        Float errornorm_ndiv = sqrt( errornorm_ndiv_sq );
                        Float errornorm_sol  = sqrt( errornorm_sol_sq  );
                        Float errornorm_curl = sqrt( errornorm_curl_sq );
                        
                        LOG << "div  error sq: " << errornorm_ndiv_sq << endl;
                        LOG << "sol  error sq: " << errornorm_sol_sq << endl;
                        LOG << "curl error sq: " << errornorm_curl_sq << endl;
                        
                        LOG << "div  error: " << errornorm_ndiv << endl;
                        LOG << "sol  error: " << errornorm_sol << endl;
                        LOG << "curl error: " << errornorm_curl << endl;
                        
                        LOG << "residual:   " << residualnorm  << endl;

                        contable << errornorm_ndiv;
                        contable << errornorm_sol;
                        contable << errornorm_curl;
                        contable << residualnorm;
                        contable << nl;

                        contable.lg();
                        


                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

                contable << nl;
                
                contable.lg();
        
            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}




