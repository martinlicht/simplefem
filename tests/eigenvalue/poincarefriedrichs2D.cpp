

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/neumannestimate.hpp"
#include "../../mesh/spanning.hpp"
#include "../../mesh/shelling.hpp"
#include "../../mesh/shelling2.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.avgsullivan.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: 2D PF estimates using shellings" << nl;
    
    LOG << "Initial mesh..." << nl;

    const bool do_gradient   = true;
    const bool do_divergence = true;
    
    // MeshSimplicial2D M = UnitSquare2D_simple();
    // MeshSimplicial2D M = UnitSquare2D_strange14();
    MeshSimplicial2D M = LShapedDomain2D();
    // MeshSimplicial2D M = SlitDomain2D();
    // MeshSimplicial2D M = SlitDomain2D_noncentered();
    // MeshSimplicial2D M = SlitDomain2D_fivetriangles();
    M.check();
    
    LOG << M.getCoordinates().text() << nl;

    LOG << M << nl;

    LOG << "Number of triangles: " << M.count_triangles() << nl;

    Float grad_estimate = NeumannEstimate( M );

    LOG << grad_estimate << nl;

    return 0;
    
    Float PF_estimate_via_shellings[2] = { notanumber, notanumber };

    //if(false)
    for( int k = 0; k < 2; k++ )
    {
        PF_estimate_via_shellings[k] = std::numeric_limits<Float>::infinity();

        auto shellings_found = generate_shellings2( M, 0 );

        LOG << shellings_found.size() << nl;

        typedef decltype(shellings_found[0]) shelling;

        std::sort( shellings_found.begin(), shellings_found.end(), 
                    [=]( const shelling& s1, const shelling& s2 ){ return s1.weight_reflection < s2.weight_reflection; } 
                    );

        for( int t = 0; t < shellings_found.size(); t++ )
        {
            const auto& shelling = shellings_found[t];
            LOG << "k=" << k << "\t" << t << "\t:\t";
            for( const auto& s : shelling ) { LOG << s << space; }
            LOG << "\t" << shelling.weight_reflection << "\t" << shelling.weight_deformation << nl;
        }

        PF_estimate_via_shellings[k] = shellings_found.front().weight_reflection;
    
    }

            
    
            
    
        
                    


    LOG << "Estimating Poincare-Friedrichs constants" << nl;

    const int min_l = 4; 
    const int max_l = 4;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
    std::vector<ConvergenceTable> contables_scalar(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables_scalar[r-min_r].table_name = "SCALAR: Mass error and numerical residuals r=" + std::to_string(r);
        contables_scalar[r-min_r] << "eigenvalue (iterated)" << "eigenvalue (post)" << "PF (post)" << "PF ratio" << "u_mass" << "du_mass" << "time" << nl;
    } 

    std::vector<ConvergenceTable> contables_vector(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables_vector[r-min_r].table_name = "VECTOR: Mass error and numerical residuals r=" + std::to_string(r);
        contables_vector[r-min_r] << "eigenvalue (iterated)" << "eigenvalue (post)" << "PF (post)" << "PF ratio" << "u_mass" << "du_mass" << "time" << nl;
    } 

    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
      
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

            LOG << "... assemble shared inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
            SparseMatrix volume_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble shared algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();


            LOG << "averaging matrix for scalar fields" << nl;

            SparseMatrix scalar_averaging = FEECSullivanAveragingMatrix( M, M.getinnerdimension(), 0, r+1, FEECAveragingMode::weighted_uniformly );


            if( do_gradient ) {
                
                LOG << "... compose system matrices (SCALAR)" << nl;

                // PING;
                // scalar_diffmatrix & scalar_incmatrix;
                // PING;
                // vector_elevationmatrix & scalar_diffmatrix;
                // PING;
                // vector_massmatrix & scalar_elevationmatrix;
                // PING;
                auto mat_A  = scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix;
                mat_A.sortandcompressentries();
                    
                LOG << "... compose CSR system matrices (SCALAR)" << nl;
        
                auto A  = MatrixCSR( mat_A  );
                
                            

                LOG << "...begin inverse iteration" << nl;
                
                const int max_attempts = 1;

                for( int s = 0; s < max_attempts; s++ )
                {

                    FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                    candidate.random(); 
                    candidate = A * candidate;
                    candidate.normalize( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    
                    const int max_inverseiterations = 10;

                    Float newratio = -1;
                    
                    FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r+1, 
                    [](const FloatVector& vec) -> FloatVector { return FloatVector({ 1. }); }
                    );

                    FloatVector conforming_one = scalar_averaging * interpol_one;

                    Float mass_of_conforming_one  = ( scalar_massmatrix * interpol_one ) * interpol_one;

                    timestamp start = timestampnow();
                    
                    for( int t = 0; t <= max_inverseiterations; t++ )
                    {
                        
                        
                        // find the next candidate

                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * candidate ) * interpol_one;
                            candidate = candidate - ( average / mass_of_conforming_one ) * conforming_one;
                        }

                        FloatVector sol( A.getdimout(), 0. );
                        
                        FloatVector rhs_sol = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate;
                        
                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * rhs_sol ) * interpol_one;
                            rhs_sol = rhs_sol - ( average / mass_of_conforming_one ) * conforming_one;
                        }

                        auto residual = sol;
                        
                        // ConjugateResidualSolverCSR( 
                        //     sol.getdimension(), 
                        //     sol.raw(), 
                        //     rhs_sol.raw(), 
                        //     A.getA(), A.getC(), A.getV(),
                        //     residual.raw(),
                        //     desired_precision,
                        //     -1
                        // );

                        {
                            DenseMatrix Bt( A.getdimout(), 1, 1. );
                            DenseMatrix B = Transpose(Bt);
                            DenseMatrix C(1,1,0.);
                            
                            auto aux     = FloatVector(1,0.);
                            auto rhs_aux = FloatVector(1,0.);
                            
                            BlockHerzogSoodhalterMethod( 
                                sol, 
                                aux, 
                                rhs_sol, 
                                rhs_aux, 
                                A, Bt, B, C, 
                                desired_precision * std::sqrt(desired_precision),
                                -1,
                                IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                            );

                        }

                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * sol ) * interpol_one;
                            sol = sol - ( average / mass_of_conforming_one ) * conforming_one;
                        }
                        
                        candidate = sol;
                        
                        
                        // assess this new candidate 

                        const auto A_candidate = A * candidate;
                        const auto M_candidate = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate; 
                        
                        const auto candidate_A_product = candidate * A_candidate; 
                        const auto candidate_M_product = candidate * M_candidate; 

                        newratio = candidate_A_product / candidate_M_product;

                        candidate /= std::sqrt(candidate_M_product); // Optional step

                        LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;
                        
                        Float u_residualmass_sq   = ( A * sol - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                        
                        LOG << "current residual: " << u_residualmass_sq << nl;

                        
                    }

                    timestamp end = timestampnow();

                    // ... computed the solution

                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    LOG << "...compute error and residual" << nl;

                    auto sol = candidate; 

                    Float u_massnorm     = sol * ( scalar_incmatrix_t * ( scalar_massmatrix * scalar_incmatrix * sol  ) );
                    Float ugrad_massnorm = sol * ( mat_A * sol );
                    Float curratio       = ugrad_massnorm / u_massnorm;
                    
                    LOG << "ratio:       " << newratio << nl;
                    LOG << "ratio:       " << curratio << nl;
                    LOG << "u mass:      " << u_massnorm << nl;
                    LOG << "u grad mass: " << ugrad_massnorm << nl;

                    LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                    
                    const Float true_eigenvalue = 1.;
                    // 1.0 is the true value 
                    // 3.0 is the true value 

                    contables_scalar[r-min_r] << newratio;
                    contables_scalar[r-min_r] << curratio;
                    contables_scalar[r-min_r] << std::sqrt(curratio);
                    contables_scalar[r-min_r] << PF_estimate_via_shellings[1] * std::sqrt(curratio);
                    contables_scalar[r-min_r] << u_massnorm;
                    contables_scalar[r-min_r] << ugrad_massnorm;
                    contables_scalar[r-min_r] << Float( end - start );
                    contables_scalar[r-min_r] << nl;

                    contables_scalar[r-min_r].lg();
                
                }
                
            } // scalar eigenvalues 








































































            if( do_divergence ){
                
                LOG << "... compose system matrices (VECTOR)" << nl;

                auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();

                PING;
                    
                auto mat_Bt = vector_incmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                PING;
                    
                auto mat_B = mat_Bt.getTranspose(); //vector_incmatrix_t & vector_massmatrix & vector_diffmatrix & vector_incmatrix; // lower left
                mat_B.sortandcompressentries();
                
                LOG << "... compose CSR system matrices (VECTOR)" << nl;
        
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                // TODO: develop preconditioners 
                // auto PA = IdentityMatrix( A.getdimin() );
                // auto PC = IdentityMatrix( C.getdimin() );

                auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                                + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix );
                auto PC = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )
                                + MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix );
                // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                            

                LOG << "...begin inverse iteration (VECTOR)" << nl;
                
                const int max_attempts = 1;

                for( int s = 0; s < max_attempts; s++ )
                {

                    FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                    candidate.random(); 
                    candidate = A * candidate;
                    candidate.normalize( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    
                    const int max_inverseiterations = 10;

                    Float newratio = -1;
                    
                    timestamp start = timestampnow();

                    for( int t = 0; t < max_inverseiterations; t++ )
                    {
                        
                        // find the next candidate

                        FloatVector sol( A.getdimout(), 0. );
                        FloatVector aux( B.getdimout(), 0. );

                        const FloatVector rhs_sol = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate;
                        const FloatVector rhs_aux = FloatVector( B.getdimout(), 0. );

                        // auto residual = sol;
                        // ConjugateResidualSolverCSR_textbook( 
                        //     sol.getdimension(), 
                        //     sol.raw(), 
                        //     rhs_sol.raw(), 
                        //     A.getA(), A.getC(), A.getV(),
                        //     residual.raw(),
                        //     desired_precision,
                        //     -1
                        // );
                        
                        // const auto PAinv = pinv(PA,desired_precision,-1);
                        // const auto PCinv = pinv(PC,desired_precision,-1);
                        // sol = PAinv * rhs_sol;
                        
                        const auto PAinv = inv(PA,desired_precision,-1);
                        const auto PCinv = inv(PC,desired_precision,-1);
                        BlockHerzogSoodhalterMethod( 
                            sol, 
                            aux, 
                            rhs_sol, 
                            rhs_aux, 
                            A, Bt, B, C, 
                            desired_precision * std::sqrt(desired_precision),
                            -1,
                            // IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                            PAinv, PCinv
                        );

                        candidate = sol;
                        
                        
                        // assess the current candidate 

                        const auto A_candidate = A * candidate;
                        const auto M_candidate = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate; 
                        
                        const auto candidate_A_product = candidate * A_candidate; 
                        const auto candidate_M_product = candidate * M_candidate; 

                        newratio = candidate_A_product / candidate_M_product;

                        candidate /= std::sqrt(candidate_M_product); // Optional step

                        LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                        Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                        Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                        
                        LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                        
                    }

                    timestamp end = timestampnow();

                    // ... computed the solution

                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    LOG << "...compute error and residual" << nl;

                    auto sol = candidate; 

                    Float u_massnorm     = sol * ( vector_incmatrix_t * ( vector_massmatrix * vector_incmatrix * sol  ) );
                    Float udiv_massnorm = sol * ( mat_A * sol );
                    Float curratio       = udiv_massnorm / u_massnorm;
                    Float u_defectmass   = ( B * sol ).norm_sq( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    
                    LOG << "ratio:           " << newratio << nl;
                    LOG << "ratio:           " << curratio << nl;
                    LOG << "u mass:          " << u_massnorm << nl;
                    LOG << "u div mass       " << udiv_massnorm << nl;
                    LOG << "u defect mass:   " << u_defectmass << nl;
                    
                    LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                    
                    const Float true_eigenvalue = 3.; // 3.0 is the true value 

                    contables_vector[r-min_r] << newratio;
                    contables_vector[r-min_r] << curratio;
                    contables_vector[r-min_r] << std::sqrt(curratio);
                    contables_vector[r-min_r] << PF_estimate_via_shellings[1] * std::sqrt(curratio);
                    contables_vector[r-min_r] << u_massnorm;
                    contables_vector[r-min_r] << udiv_massnorm;
                    contables_vector[r-min_r] << Float( end - start );
                    contables_vector[r-min_r] << nl;

                    contables_vector[r-min_r].lg();
                
                }
                
            } // vector eigenvalue 

            
        } // poly degree

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        

    } // levels

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}


