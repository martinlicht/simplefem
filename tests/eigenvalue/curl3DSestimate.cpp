

/**/

#include <cmath>
#include <cstdlib>

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
#include "../../fem/global.inclsullivan.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D curl estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    
    
    
    
    
    // ================================================================================
    // We have different choices of meshes 
    // ================================================================================
    
    std::vector<MeshSimplicial3D> Ms = 
    {
        UnitSimplex3D(),    // 0 Unit simplex Dirichlet from 0
        UnitSimplex3D(),    // 1 Unit simplex Dirichlet from 1 
        UnitSimplex3D(),    // 2 Unit simplex Dirichlet from 2
        UnitSimplex3D(),    // 3 Unit simplex Dirichlet from 3
        UnitSimplex3D(),    // 4 Unit simplex Neumann
        UnitSimplex3D(),    // 5 Distorted simplex 
        UnitCube3D(),       // 6 Unit cube 
        UnitCube3D(),       // 7 Unit cube (Dirichlet)
        FicheraCorner3D(),  // 8 Fichera corner
        CrossedBricks3D()   // 9 Crossed bricks 
    };

    for( int f = 0; f <= 3; f++ ) { Ms[0].set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet ); }
    for( int f = 1; f <= 3; f++ ) { Ms[1].set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet ); }
    for( int f = 2; f <= 3; f++ ) { Ms[2].set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet ); }
    for( int f = 3; f <= 3; f++ ) { Ms[3].set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet ); }
    // Nothing for #4, which has Neumann boundary conditions
    for( int i = 0; i <= 4; i++ ) { Ms[i].complete_dirichlet_flags_from_facets(); Ms[i].check_dirichlet_flags(false); }

    Ms[5].getCoordinates().scale( 3. );
    Ms[5].getCoordinates().setdata(3, 2, 1./3. );

    Ms[7].automatic_dirichlet_flags();
    
    struct ReferencePair { Float eigenvalue; Float PF_estimate; };

    // TODO(martin): fill in the correct values here 
    
    const Float PF_factor = 2. / Constants::pi;

    std::vector<ReferencePair> curl_reference_values = 
    {
        { notanumber, PF_factor * Ms[0].getMeshDiameter() },    // 0 Unit simplex Dirichlet from 0
        { notanumber, PF_factor * Ms[1].getMeshDiameter() },    // 1 Unit simplex Dirichlet from 1 
        { notanumber, PF_factor * Ms[2].getMeshDiameter() },    // 2 Unit simplex Dirichlet from 2
        { notanumber, PF_factor * Ms[3].getMeshDiameter() },    // 3 Unit simplex Dirichlet from 3
        { notanumber, PF_factor * Ms[4].getMeshDiameter() },    // 4 Unit simplex Neumann
        { notanumber, PF_factor * Ms[5].getMeshDiameter() },    // 5 Distorted simplex 
        { notanumber, PF_factor * Ms[6].getMeshDiameter() },    // 6 Unit cube 
        { notanumber, PF_factor * Ms[7].getMeshDiameter() },    // 7 Unit cube (Dirichlet)
        { notanumber,                          notanumber },    // 8 Fichera corner
        { notanumber,                          notanumber }     // 9 Crossed bricks 
    };

    // ================================================================================
    // Select the mesh based on the input 
    // ================================================================================
    
    assert( Ms.size() == curl_reference_values.size() );
    
    unsigned int choice_of_mesh = 0;
    
    if( argc > 1 )
    {
        const char* end = nullptr;
        bool has_overflown;
        int value = string_to_integer( argv[1], &end, 10, has_overflown );
        assert( end != nullptr );

        // Check if the entire argument was parsed and within int range
        if( *end != '\0' ) {

            LOG << "Error: The provided argument is not a valid integer:" << argv[1] << "\n";

        } else if( value != static_cast<int>(value) ) {

            LOG << "Error: The provided argument is out of 'int' range.\n";

        } else if( value < 0 || value >= Ms.size() ) {

            LOG << "Error: The provided argument is not a mesh index: " << value << "\n";

        } else {
            
            choice_of_mesh = static_cast<int>(value);
            LOG << "Select mesh based on input: " << choice_of_mesh << nl;

        }
    } else {
        LOG << "Default mesh selected: " << choice_of_mesh << nl;
    }

    assert( choice_of_mesh < Ms.size() );

    // ================================================================================
    // Finally, we select a mesh
    // ================================================================================
    
    MeshSimplicial3D& M = Ms[ choice_of_mesh ]; 
    
    M.check();
    M.check_dirichlet_flags(false);

    
    
    
    
    

    // ================================================================================
    // Main loop
    // ================================================================================
    
    LOG << "Estimating Poincare-Friedrichs constant of curl operator (Sullivan)" << nl;

    const int min_l = 2; 
    const int max_l = 2;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
    std::vector<ConvergenceTable> contables(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables[r-min_r] << "PF" << "computed" << "ratioE" << "log_2(ratioE)" << "ratioP" << "log_2(ratioP)" << "u_mass" << "du_mass" << "time" << nl;
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
    
            auto scalar_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 ) );
            auto vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   ) );
            auto pseudo_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ) );

            LOG << "... assemble inclusion matrices" << nl;
    
            auto scalar_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 ) );
            auto scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            auto vector_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   ) );
            auto vector_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            auto scalar_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 ) );
            auto scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            auto vector_diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r ) );
            auto vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = Conjugation( pseudo_massmatrix, vector_diffmatrix & vector_incmatrix );
            // mat_A.sortandcompressentries();
                
            auto mat_Bt = ( vector_incmatrix_t & vector_massmatrix ) & ( scalar_diffmatrix & scalar_incmatrix ); // upper right
            // mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & pseudo_massmatrix & diffmatrix & vector_incmatrix; // lower left
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
                      Conjugation( vector_massmatrix, scalar_diffmatrix & scalar_incmatrix );
            
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
                
                LOG << "ratio:           " << newratio << nl;
                LOG << "ratio:           " << curratio << nl;
                LOG << "u mass:          " << u_massnorm << nl;
                LOG << "u curl mass      " << ucurl_massnorm << nl;
                LOG << "u defect mass:   " << u_defectmass << nl;
                
                LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                
                // const Float true_eigenvalue = 2.; // 3.0 is the true value 

                const Float reference_eigenvalue = curl_reference_values[choice_of_mesh].eigenvalue;

                const Float reference_pf = curl_reference_values[choice_of_mesh].PF_estimate;

                // contables[r-min_r] << "PF" << "computed" << "ratioE" << "log_2(ratioE)" << "ratioP" << "log_2(ratioP)" << "u_mass" << "du_mass" << "time" << nl;
                
                contables[r-min_r] << std::sqrt(newratio);
                contables[r-min_r] << newratio;
                contables[r-min_r] << newratio / reference_eigenvalue;
                contables[r-min_r] << std::log2( newratio / reference_eigenvalue ); 
                contables[r-min_r] << std::sqrt(newratio) / reference_pf;
                contables[r-min_r] << std::log2( std::sqrt(newratio) / reference_pf ); 
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
