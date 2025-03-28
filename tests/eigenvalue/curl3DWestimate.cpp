

/**/

#include <cmath>

#include <string>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.inclwhitney.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D curl estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitCube3D();
    M.getCoordinates().scale( Constants::pi );
    
    M.check();
    
    // M.automatic_dirichlet_flags();
    // M.check_dirichlet_flags();

    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    
            
    
    
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of curl operator (Whitney)" << nl;

    const int min_l = 0; 
    const int max_l = 3;
    
    const int min_r = 1;
    const int max_r = 2;
    
    
    std::vector<ConvergenceTable> contables(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables[r-min_r] << "eigenvalue" << "ratio" << "log_2(ratio)" << "diff" << "log_2(diff)" << "u_mass" << "du_mass" << "time" << nl;
    } 

    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
      
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Level: "             << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "Polynomial degree: " << min_r << " <= " << r << " <= " << max_r << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            auto scalar_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   ) );
            auto vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   ) );
            auto pseudo_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ) );

            LOG << "... assemble inclusion matrices" << nl;
    
            auto scalar_incmatrix   = MatrixCSR( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r ) );
            auto scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            auto vector_incmatrix   = MatrixCSR( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r ) );
            auto vector_incmatrix_t = vector_incmatrix.getTranspose();


            LOG << "... assemble algebraic matrices" << nl;
    
            auto scalar_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r ) );
            auto scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            auto vector_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r ) );
            auto vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            auto vector_elevationmatrix = MatrixCSR( FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 ) );
            auto vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                    

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = Conjugation( pseudo_massmatrix, vector_diffmatrix & vector_incmatrix );
            // mat_A.sortandcompressentries();
                
            auto mat_Bt = vector_incmatrix_t & vector_massmatrix & ( vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix ); // upper right
            // mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix; // lower left
            // mat_B.sortandcompressentries();
            
            LOG << "... compose CSR system matrices" << nl;
    
            // auto A  = MatrixCSR( mat_A  );
            // auto Bt = MatrixCSR( mat_Bt );
            // auto B  = MatrixCSR( mat_B  );
            auto& A  = mat_A;
            auto& Bt = mat_Bt;
            auto& B  = mat_B;
            
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );

            auto PA = Conjugation( vector_massmatrix, vector_incmatrix )
                      + 
                      Conjugation( pseudo_massmatrix, vector_diffmatrix & vector_incmatrix );
            
            auto PC = Conjugation( scalar_massmatrix, scalar_incmatrix )
                      + 
                      Conjugation( vector_massmatrix, vector_elevationmatrix & ( scalar_diffmatrix & scalar_incmatrix ) );
            
            // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        

            LOG << "... begin inverse iteration" << nl;
            
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

                    Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    
                    LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                    
                }

                timestamp end = timestampnow();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "... compute error and residual" << nl;

                auto sol = candidate; 

                Float u_massnorm     = sol * ( vector_incmatrix_t * ( vector_massmatrix * vector_incmatrix * sol  ) );
                Float ucurl_massnorm = sol * ( mat_A * sol );
                Float curratio       = ucurl_massnorm / u_massnorm;
                Float u_defectmass   = ( B * sol ).norm_sq( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                
                LOG << "ratio:         " << newratio << nl;
                LOG << "ratio:         " << curratio << nl;
                LOG << "u mass:        " << u_massnorm << nl;
                LOG << "u curl mass    " << ucurl_massnorm << nl;
                LOG << "u defect mass: " << u_defectmass << nl;
                
                const Float true_eigenvalue = 2.; // 3.0 is the true value 

                contables[r-min_r] << newratio;
                contables[r-min_r] << newratio / true_eigenvalue - 1.;
                contables[r-min_r] << - std::log2( newratio / true_eigenvalue - 1. ); 
                contables[r-min_r] << newratio - true_eigenvalue;
                contables[r-min_r] << std::log2( newratio - true_eigenvalue );
                contables[r-min_r] << u_massnorm;
                contables[r-min_r] << ucurl_massnorm;
                contables[r-min_r] << Float( end - start );
                contables[r-min_r] << nl;

                contables[r-min_r].lg();
            
            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

    } 

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
