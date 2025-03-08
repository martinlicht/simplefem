

/**/

#include <cmath>
#include <cstdlib>

#include <string>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: 3D divergence estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitSimplex3D(); 
    // MeshSimplicial3D M = UnitCube3D(); M.getCoordinates().scale( Constants::pi );
    
    // TODO: This is to get a badly-shaped simplex 
    M.getCoordinates().scale( 3. ); M.getCoordinates().setdata(3, 2, 1./3. );
    
    M.check();
    
    // M.automatic_dirichlet_flags();
    // M.check_dirichlet_flags();

    {
        LOG << "Fine-tuned boundary conditions" << nl;

        int first_bc_face = 0;
        
        if( argc > 1 )
        {
            // Parse using strtol
            char* end = nullptr;
            long val = std::strtol(argv[1], &end, 10);

            // Check if the entire argument was parsed and within int range
            if( *end != '\0' ) {

                LOG << "Error: The provided argument is not a valid integer:" << argv[1] << "\n";

            } else if( val != static_cast<int>(val) ) {

                LOG << "Error: The provided argument is out of 'int' range.\n";

            } else {
                
                int number = static_cast<int>(val);
                LOG << "Command-line argument provided: " << number << nl;

                first_bc_face = number;

            }
        } else {
            LOG << "Dirichlet faces start, per default, at " << first_bc_face << nl;
        }

        for( int f = first_bc_face; f <= 3; f++ ) {
            M.set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet );
        }
        M.complete_dirichlet_flags_from_facets();
        M.check_dirichlet_flags(false);
    
    }

    for( int f = 0; f < M.count_faces(); f++ )
    {
        LOG << f << ": ";
        for( int i = 0; i <= 2; i++ ) LOG << M.get_subsimplex( 2, 0, f, i ) << space;
        LOG << " > " << (int)M.get_flag( 2, f ) << nl;
    }
    
    
    
            
    
    
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of div operator (Whitney)" << nl;

    const int min_l = 4; 
    const int max_l = 4;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
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
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r   );
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 );
            
            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix pseudo_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r   );
            SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            SparseMatrix pseudo_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r );
            SparseMatrix pseudo_diffmatrix_t = pseudo_diffmatrix.getTranspose();

            SparseMatrix pseudo_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
            SparseMatrix pseudo_elevationmatrix_t = pseudo_elevationmatrix.getTranspose();


            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = pseudo_incmatrix_t & pseudo_diffmatrix_t & volume_massmatrix & pseudo_diffmatrix & pseudo_incmatrix;
            mat_A.sortandcompressentries();
                
            auto mat_Bt = pseudo_incmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix; // upper right
            mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //pseudo_incmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix; // lower left
            mat_B.sortandcompressentries();
            
            LOG << "... compose CSR system matrices" << nl;
    
            auto A  = MatrixCSR( mat_A  );
            auto Bt = MatrixCSR( mat_Bt );
            auto B  = MatrixCSR( mat_B  );
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // TODO: develop preconditioners 
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );

            auto PA = MatrixCSR( pseudo_incmatrix_t & pseudo_massmatrix & pseudo_incmatrix )
                      +
                      MatrixCSR( pseudo_incmatrix_t & pseudo_diffmatrix_t & volume_massmatrix & pseudo_diffmatrix & pseudo_incmatrix );
            auto PC = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                      +
                      MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & pseudo_elevationmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix );
            // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        

            LOG << "...begin inverse iteration" << nl;
            
            const int max_attempts = 1;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                candidate.random(); 
                candidate = A * candidate;
                candidate.normalize( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ); 
                
                const int max_inverseiterations = 10;

                Float newratio = -1;
                
                timestamp start = timestampnow();

                for( int t = 0; t < max_inverseiterations; t++ )
                {
                    
                    
                    // find the next candidate

                    FloatVector sol( A.getdimout(), 0. );
                    FloatVector aux( B.getdimout(), 0. );

                    const FloatVector rhs_sol = ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ) * candidate;
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
                    const auto M_candidate = ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ) * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    candidate /= std::sqrt(candidate_M_product); // Optional step

                    LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ); 
                    
                    LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                    
                }

                timestamp end = timestampnow();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual" << nl;

                auto sol = candidate; 

                Float u_massnorm     = sol * ( pseudo_incmatrix_t * ( pseudo_massmatrix * pseudo_incmatrix * sol  ) );
                Float udiv_massnorm = sol * ( mat_A * sol );
                Float curratio       = udiv_massnorm / u_massnorm;
                Float u_defectmass   = ( B * sol ).norm_sq( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                
                LOG << "ratio:         " << newratio << nl;
                LOG << "ratio:         " << curratio << nl;
                LOG << "u mass:        " << u_massnorm << nl;
                LOG << "u div mass     " << udiv_massnorm << nl;
                LOG << "u defect mass: " << u_defectmass << nl;
                
                LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                
                const Float true_eigenvalue = 3.; // 3.0 is the true value 

                contables[r-min_r] << newratio;
                contables[r-min_r] << newratio / true_eigenvalue - 1.;
                contables[r-min_r] << - std::log2( newratio / true_eigenvalue - 1. ); 
                contables[r-min_r] << newratio - true_eigenvalue;
                contables[r-min_r] << std::log2( newratio - true_eigenvalue );
                contables[r-min_r] << u_massnorm;
                contables[r-min_r] << udiv_massnorm;
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
