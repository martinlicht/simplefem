

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
    
    contable << "ratio" << "u_mass" << "du_mass" << "time";
    

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
            
            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = scalar_incmatrix_t & scalar_diffmatrix_t & ( vector_massmatrix ) & scalar_diffmatrix & scalar_incmatrix;
            mat_A.sortandcompressentries();
                
            LOG << "... compose CSR system matrices" << nl;
    
            auto A  = MatrixCSR( mat_A  );
            
                        

            LOG << "...begin inverse iteration" << nl;
            
            const int max_attempts = 5;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                candidate.random(); 
                candidate.normalize( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                
                const int max_inverseiterations = 10;

                Float newratio = -1;

                
                
                timestamp start = gettimestamp();
                for( int t = 0; t <= max_inverseiterations; t++ )
                {
                    
                    
                    // assess the current candidate 

                    const auto A_candidate = A * candidate;
                    const auto M_candidate = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    LOG << "current ratio : " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    
                    // find the next candidate

                    FloatVector sol( A.getdimout(), 0. );
                    
                    const FloatVector rhs_sol = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate / sqrt( candidate_M_product );
                    
                    auto residual = sol;
                    ConjugateGradientSolverCSR( 
                        sol.getdimension(), 
                        sol.raw(), 
                        rhs_sol.raw(), 
                        A.getA(), A.getC(), A.getV(),
                        residual.raw(),
                        desired_precision / 10000,
                        -1
                    );
                    
                    
                    candidate = sol;
                    
                    
                    Float u_residualmass_sq   = ( A * sol - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    LOG << "current residual : " << u_residualmass_sq << nl;

                    
                }

                timestamp end = gettimestamp();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual:" << nl;

                auto sol = candidate; 

                Float u_massnorm     = sol * ( scalar_incmatrix_t * ( scalar_massmatrix * scalar_incmatrix * sol  ) );
                Float ugrad_massnorm = sol * ( mat_A * sol );
                Float curratio       = ugrad_massnorm / u_massnorm;
                
                LOG << "ratio:           " << newratio << nl;
                LOG << "ratio:           " << curratio << nl;
                LOG << "u mass:          " << u_massnorm << nl;
                LOG << "u grad mass      " << ugrad_massnorm << nl;
                
                contable << newratio;
                contable << u_massnorm;
                contable << ugrad_massnorm;
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
