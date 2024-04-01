

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
#include "../../fem/global.cechmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: Cech complex 3D estimates" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitSimplex3D();
    // MeshSimplicial3D M = UnitCube3D();
    
    M.check();
    
    // M.automatic_dirichlet_flags();
    // M.check_dirichlet_flags();

    LOG << "Estimating Poincare-Friedrichs constant of Cech complex" << nl;

    const int min_l = 0; 
    const int max_l = 6;

    const int n = M.getinnerdimension();
    
    ConvergenceTable contable;
    contable.table_name = "Cech complex estimates";
    contable << "level" << "#V" << "#E" << "#F" << "#T" << "C_grad" << "C_curl" << "C_div" << nl;
    
    assert( 0 <= min_l and min_l <= max_l );
    
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;

        contable << (float)l;
        contable << (float)M.count_vertices() << (float)M.count_edges() << (float)M.count_faces() << (float)M.count_tetrahedra();

        LOG << "... assemble matrices" << nl;

        SparseMatrix vertexones = SparseMatrix( M.count_vertices(), 1, M.count_vertices(), [](int r)->SparseMatrix::MatrixEntry{ return SparseMatrix::MatrixEntry(r,0,1.0); } );
    
        std::vector<SparseMatrix> cech_massmatrix;
        for( int k = 0; k <= n; k++ ) cech_massmatrix.push_back( FEECCechMassMatrix( M, M.getinnerdimension(), k, -k ) );

        std::vector<SparseMatrix> cech_diffmatrix;
        for( int k = 0; k <= n-1; k++ ) cech_diffmatrix.push_back( FEECCechDiffMatrix( M, M.getinnerdimension(), k ) );

        std::vector<SparseMatrix> cech_diffmatrix_t;
        for( int k = 0; k <= n-1; k++ ) cech_diffmatrix_t.push_back( cech_diffmatrix[k].getTranspose() );

        std::vector<SparseMatrix> A;
        for( int k = 0; k <= n-1; k++ ) A.push_back( cech_diffmatrix_t[k] & cech_massmatrix[k+1] & cech_diffmatrix[k] );

        std::vector<SparseMatrix> B;
        B.push_back( vertexones.getTranspose() & cech_massmatrix[0] );
        for( int k = 1; k <= n-1; k++ ) B.push_back( cech_diffmatrix_t[k-1] & cech_massmatrix[k] );

        std::vector<SparseMatrix> Bt;
        Bt.push_back( cech_massmatrix[0] & vertexones );
        for( int k = 1; k <= n-1; k++ ) Bt.push_back( cech_massmatrix[k] & cech_diffmatrix[k-1] );

        std::vector<ZeroOperator> C;
        for( int k = 0; k <= n-1; k++ ) C.push_back( ZeroOperator( Bt[k].getdimin() ) );

        LOG << "... compute lowest eigenvalues" << nl;


        // LOG << M << nl;
        // for( const auto& mat : cech_massmatrix ) LOG << mat << nl;
        
        for( int k = 0; k <= n-1; k++ )
        {
            
            LOG << "...begin inverse iteration: " << k << nl;
            
            const int max_attempts = 1;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A[k].getdimout(), 0. ); 
                candidate.random(); 
                candidate = A[k] * candidate;
                candidate.normalize( cech_massmatrix[k] ); 
                
                const int max_inverseiterations = 3;

                Float newratio = -1;
                
                timestamp start = timestampnow();

                if( k == 0 )
                for( int t = 0; t < max_inverseiterations; t++ )
                {

                    // find the next candidate

                    FloatVector sol( A[k].getdimout(), 0. );
                    FloatVector aux( B[k].getdimout(), 0. );

                    const FloatVector rhs_sol = cech_massmatrix[k] * candidate;
                    const FloatVector rhs_aux = FloatVector( B[k].getdimout(), 0. );

                    BlockHerzogSoodhalterMethod( 
                        sol, 
                        aux, 
                        rhs_sol, 
                        rhs_aux, 
                        A[k], Bt[k], B[k], C[k], 
                        desired_precision * sqrt(desired_precision),
                        -1,
                        IdentityMatrix( A[k].getdimin() ), IdentityMatrix( C[k].getdimin() ) 
                    );

                    candidate = sol;
                    
                    
                    // assess the current candidate 

                    const auto A_candidate = A[k]               * candidate;
                    const auto M_candidate = cech_massmatrix[k] * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    candidate /= sqrt(candidate_M_product); // Optional step

                    LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    Float u_residualmass_sq   = ( A[k] * sol + Bt[k] * aux - rhs_sol ).norm_sq(); 
                    Float aux_residualmass_sq = ( B[k] * sol               - rhs_aux ).norm_sq(); 
                    
                    LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                    
                }

                timestamp end = timestampnow();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual" << nl;

                auto sol = candidate; 

                // assess the current candidate 

                const auto A_sol = A[k]               * sol;
                const auto M_sol = cech_massmatrix[k] * sol; 
                
                const auto sol_A_product = sol * A_sol; 
                const auto sol_M_product = sol * M_sol; 

                newratio = sol_A_product / sol_M_product;
            
                const Float defect_euclnorm = ( B[k] * sol ).norm_sq();
            
                contable << newratio;

                LOG << "ratio: " << newratio << " defect: " << defect_euclnorm << nl;

                
            }
            
        }

        contable << nl;

        contable.lg();

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        

    } 

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
