

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Solution of Darcy Problem" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitCube3D();
    M.getcoordinates().scale( Constants::pi );
    
    M.check();
    
    M.automatic_dirichlet_flags();
    M.check_dirichlet_flags();

    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    const Float A = Constants::twopi;

            
    
    
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of curl operator" << nl;

    const int min_l = 0; 
    const int max_l = 7;
    
    const int min_r = 1;
    const int max_r = 1;

    const int r_plus_vector = 0;
    const int r_plus_scalar = 0;
    const int r_plus_pseudo = 0;
    
    
    ConvergenceTable contable("Mass error and numerical residuals");
    
    contable << "ratio" << "u_mass" << "du_mass" << "sigma_mass" << "u_defect" << "aux_defect" << "time";
    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    assert( 0 <= r_plus_scalar and 0 <= r_plus_vector and 0 <= r_plus_pseudo );
      
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "integration with: " << r_plus_vector << ", " << r_plus_scalar << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 + r_plus_scalar );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   + r_plus_vector );
            SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 + r_plus_pseudo );

            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix;
            mat_A.sortandcompressentries();
                
            auto mat_Bt = vector_incmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
            mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & pseudo_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
            mat_B.sortandcompressentries();
            
            LOG << "... compose CSR system matrices" << nl;
    
            auto A  = MatrixCSR( mat_A  );
            auto Bt = MatrixCSR( mat_Bt );
            auto B  = MatrixCSR( mat_B  );
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // TODO: develop preconditioners 
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );

            auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                              + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix );
            auto PC = MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix );
            // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        

            LOG << "...begin inverse iteration" << nl;
            
            const int max_attempts = 5;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                candidate.random(); 
                candidate.normalize( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                
                const int max_inverseiterations = 10;

                Float newratio = -1;

                
                
                const FloatVector rhs_aux = FloatVector( B.getdimout(), 0. );
                
                // FloatVector sol( A.getdimout(), 0. );
                FloatVector aux( B.getdimout(), 0. );

                timestamp start = gettimestamp();
                for( int t = 0; t < max_inverseiterations; t++ )
                {
                    
                    
                    // assess the current candidate 

                    const auto A_candidate = A * candidate;
                    const auto M_candidate = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    LOG << "current ratio : " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    
                    // find the next candidate

                    FloatVector sol( A.getdimout(), 0. );
                    FloatVector aux( B.getdimout(), 0. );

                    const FloatVector rhs_sol = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate / sqrt( candidate_M_product );
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
                    
                    // const auto PAinv = pinv(PA,desired_precision,-1);
                    // const auto PCinv = pinv(PC,desired_precision,-1);
                    BlockHerzogSoodhalterMethod( 
                        sol, 
                        aux, 
                        rhs_sol, 
                        rhs_aux, 
                        A, Bt, B, C, 
                        desired_precision * sqrt(desired_precision),
                        -1,
                        IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                        // PAinv, PCinv
                    );

                    candidate = sol;
                    
                    // compute the relevant ratio and next right-hand side
                    
                    // rhs_sol = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * sol;
                    
                    // Float sol_massnorm_sq = rhs_sol * sol;

                    // rhs_sol /= sqrt( sol_massnorm_sq );

                    // Float sol_Anorm_sq = sol * ( A * sol );
                    
                    // newratio = sqrt( sol_Anorm_sq / sol_massnorm_sq );
                    
                    Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    LOG << "current residuals : " << u_residualmass_sq << " " << aux_residualmass_sq << nl;

                    
                }

                timestamp end = gettimestamp();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual:" << nl;

                auto sol = candidate; 

                Float u_massnorm     = sqrt( sol * ( vector_incmatrix_t * ( vector_massmatrix * vector_incmatrix * sol  ) ) );
                Float ucurl_massnorm = sqrt( sol * ( mat_A * sol ) );
                Float aux_massnorm   = sqrt( aux * ( scalar_incmatrix_t * ( scalar_massmatrix * scalar_incmatrix * aux  ) ) );
                Float u_defectmass   = ( B * sol ).norm( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                Float aux_defectmass = ( Bt * aux ).norm( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                
                LOG << "ratio:           " << newratio << nl;
                LOG << "u mass:          " << u_massnorm << nl;
                LOG << "u curl mass      " << ucurl_massnorm << nl;
                LOG << "aux mass:        " << aux_massnorm << nl;
                LOG << "u defect mass:   " << u_defectmass << nl;
                LOG << "aux defect mass: " << aux_defectmass << nl;
                
                contable << newratio;
                contable << u_massnorm;
                contable << ucurl_massnorm;
                contable << aux_massnorm;
                contable << u_defectmass;
                contable << aux_defectmass;
                contable << Float( end - start );
                contable << nl;

                contable.lg();
            
            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        
        

    } 

    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
