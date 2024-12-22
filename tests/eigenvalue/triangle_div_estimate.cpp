

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
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

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: 2D grad estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D M = UnitTriangle2D();
    M.getCoordinates().scale( Constants::pi );
    
    M.check();
    
    // M.automatic_dirichlet_flags();
    auto f0 = SimplexFlag::SimplexFlagDirichlet;
    auto f1 = SimplexFlag::SimplexFlagDirichlet;
    auto f2 = SimplexFlag::SimplexFlagDirichlet;
    M.set_flags( 1, { f0, f1, f2 } );

    M.check_dirichlet_flags();

    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    
            
    
    
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of grad operator (Sullivan)" << nl;

    const int min_l = 3; 
    const int max_l = 3;
    
    const int min_r = 2;
    const int max_r = 2;
    
    
    std::vector<ConvergenceTable> contables(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables[r-min_r] << "eigenvalue" << "u_mass" << "du_mass" << "time" << nl;
    } 
    
    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
      
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            
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
            
            const int max_attempts = 1;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                candidate.random(); 
                candidate = A * candidate;
                LOG << candidate.norm() << space << candidate.getdimension() << nl;
                candidate.normalize( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                
                const int max_inverseiterations = 10;

                Float newratio = -1;
                
                timestamp start = timestampnow();
                
                for( int t = 0; t <= max_inverseiterations; t++ )
                {
                    
                    
                    // find the next candidate

                    FloatVector sol( A.getdimout(), 0. );
                    
                    FloatVector rhs_sol = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate;
                    
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
                            desired_precision * sqrt(desired_precision),
                            -1,
                            IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                        );

                    }

                    candidate = sol;
                    
                    
                    // assess this new candidate 

                    const auto A_candidate = A * candidate;
                    const auto M_candidate = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    candidate /= sqrt(candidate_M_product); // Optional step

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
                
                contables[r-min_r] << newratio;
                contables[r-min_r] << u_massnorm;
                contables[r-min_r] << ugrad_massnorm;
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
