

/**/

#include <cmath>

#include <string>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.cechmatrix.hpp"


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
    const int max_l = 4;

    const int n = M.getinnerdimension();
    
    ConvergenceTable contable_results;
    contable_results.table_name = "Cech complex estimates";
    contable_results << "level" << "#V" << "#E" << "#F" << "#T" << "C_grad" << "C_curl" << "C_div" << nl;
    
    ConvergenceTable contable_parameters;
    contable_parameters.table_name = "Mesh parameters";
    contable_parameters << "level" << "maxD" << "shape" << "comp" << "radii" << "height" << "size" << nl;
    
    assert( 0 <= min_l and min_l <= max_l );
    
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;

        contable_results << (Float)l;
        contable_results << (Float)M.count_vertices() << (Float)M.count_edges() << (Float)M.count_faces() << (Float)M.count_tetrahedra();

        LOG << "... assemble matrices" << nl;

        SparseMatrix vertexones = SparseMatrix( M.count_vertices(), 1, M.count_vertices(), [](int r)->SparseMatrix::MatrixEntry{ return SparseMatrix::MatrixEntry(r,0,1.0); } );
    
        // std::vector<SparseMatrix> bar;
        // bar.push_back( FEECCechMassMatrix( M, M.getinnerdimension(), 1, 0 ) );

        // return 0;

        std::vector<SparseMatrix> cech_massmatrix; cech_massmatrix.reserve(4); // TODO: wierd behavior, calls copy constructor
        for( int k = 0; k <= n; k++ ) cech_massmatrix.push_back( FEECCechMassMatrix( M, M.getinnerdimension(), k, 0 ) );

        std::vector<SparseMatrix> cech_diffmatrix; cech_diffmatrix.reserve(3);
        for( int k = 0; k <= n-1; k++ ) cech_diffmatrix.push_back( FEECCechDiffMatrix( M, M.getinnerdimension(), k ) );

        std::vector<SparseMatrix> cech_diffmatrix_t; cech_diffmatrix_t.reserve(3);
        for( int k = 0; k <= n-1; k++ ) cech_diffmatrix_t.push_back( cech_diffmatrix[k].getTranspose() );

        LOG << "... assemble composed matrices" << nl;

        std::vector<SparseMatrix> A; A.reserve(3);
        for( int k = 0; k <= n-1; k++ ) A.push_back( cech_diffmatrix_t[k] & cech_massmatrix[k+1] & cech_diffmatrix[k] );

        std::vector<SparseMatrix> B; B.reserve(3);
        B.push_back( vertexones.getTranspose() & cech_massmatrix[0] );
        for( int k = 1; k <= n-1; k++ ) B.push_back( cech_diffmatrix_t[k-1] & cech_massmatrix[k] );

        std::vector<SparseMatrix> Bt; Bt.reserve(3);
        Bt.push_back( cech_massmatrix[0] & vertexones );
        for( int k = 1; k <= n-1; k++ ) Bt.push_back( cech_massmatrix[k] & cech_diffmatrix[k-1] );

        std::vector<ZeroOperator> C; C.reserve(3);
        for( int k = 0; k <= n-1; k++ ) C.push_back( ZeroOperator( Bt[k].getdimin() ) );

        
        
        
        LOG << "... compute mesh parameters" << nl;

        const Float maxdiameter         = M.getMaximumDiameter();

        const Float shape_measure       = M.getShapemeasure();

        const Float comparison_quotient = M.getComparisonQuotient();

        const Float height_quotient     = M.getHeightQuotient();
        
        const Float radii_quotient[4]   = { M.getRadiiQuotient(0), M.getRadiiQuotient(1), M.getRadiiQuotient(2), M.getRadiiQuotient(3) };
        
        const int   patch_size          = M.getVertexPatchSize(); 

        const int   supsimplex_size[4]  = { M.getSupersimplexSize(0), M.getSupersimplexSize(1), M.getSupersimplexSize(2), -1 };
        
        DenseMatrix Antihor( n+1, n+1, notanumber );
        DenseMatrix Diffver( n+1, n+1, notanumber );

        for( int i = 0; i <= n; i++ )
        for( int j = 0; j <= n; j++ )
        {
            if( j != 0 )
                // Antihor(i,j) = std::sqrt( n / ( ( 2. + n ) * j*j ) * power_numerical( radii_quotient[i] + 1., n ) ) * ( 1. + 1./radii_quotient[i] ) ;
                Antihor(i,j) = 2. / ( Constants::pi * std::sqrt((Float)j) );
            
            Antihor(0,1) = 1. / Constants::pi;
            
            if( i != 3 ) 
                // Diffver(i,j) = std::sqrt( n * supsimplex_size[i] );
                Diffver(i,j) = std::sqrt( n * (n-i) ); // TODO: full proof in the manuscript
        }


        {
            LOG << "patch size:          " << M.getVertexPatchSize() << nl;
            LOG << "maximum diameter:    " << M.getMaximumDiameter() << nl;
            LOG << "minumum diameter:    " << M.getMinimumDiameter() << nl;
            LOG << "comparison quotient: " << M.getComparisonQuotient() << nl;
            LOG << "radii quotient:      " << M.getRadiiQuotient(0) << nl;
            LOG << "height ratio:        " << M.getHeightQuotient() << nl;
            LOG << "shape measure:       " << M.getShapemeasure() << nl;
        }

        contable_parameters << (Float)l << maxdiameter << shape_measure << comparison_quotient << radii_quotient[3] << height_quotient << (Float)patch_size << nl;
        
        LOG << "... compute lowest eigenvalues" << nl;
        
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
                
                const int max_inverseiterations = 5;

                timestamp start = timestampnow();

                for( int t = 0; t < max_inverseiterations; t++ )
                {

                    LOG << "...purification: " << t << nl;
            
                    FloatVector sol( A[k].getdimout(), 0. ); sol = candidate;
                    FloatVector aux( B[k].getdimout(), 0. );

                    const FloatVector rhs_sol = cech_massmatrix[k] * candidate;
                    const FloatVector rhs_aux = FloatVector( B[k].getdimout(), 0. );

                    BlockHerzogSoodhalterMethod( 
                        sol, 
                        aux, 
                        rhs_sol, 
                        rhs_aux, 
                        A[k], Bt[k], B[k], C[k], 
                        desired_precision * std::sqrt(desired_precision),
                        0,
                        IdentityMatrix( A[k].getdimin() ), IdentityMatrix( C[k].getdimin() ) 
                    );

                    candidate = sol;
                    
                    
                    // assess the current candidate 

                    const auto A_candidate = A[k]               * candidate;
                    const auto M_candidate = cech_massmatrix[k] * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    Float rayleigh_quotient = -1;
                
                    rayleigh_quotient = candidate_A_product / candidate_M_product;

                    candidate /= std::sqrt(candidate_M_product); // Optional step

                    LOG << "current ratio: " << rayleigh_quotient << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    Float u_residualmass_sq   = ( A[k] * sol + Bt[k] * aux - rhs_sol ).norm_sq(); 
                    Float aux_residualmass_sq = ( B[k] * sol               - rhs_aux ).norm_sq(); 
                    
                    LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                    
                }

                timestamp end = timestampnow();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual" << nl;

                auto eigenvector = candidate; 

                // assess the current candidate 

                const auto A_eigenvector = A[k]               * eigenvector;
                const auto M_eigenvector = cech_massmatrix[k] * eigenvector; 
                
                const auto eigenvector_A_product = eigenvector * A_eigenvector; 
                const auto eigenvector_M_product = eigenvector * M_eigenvector; 

                Float eigenvalue = eigenvector_A_product / eigenvector_M_product;
            
                const Float defect_euclnorm = ( B[k] * eigenvector ).norm_sq();
            
                contable_results << eigenvalue;

                LOG << "ratio: " << eigenvalue << " defect: " << defect_euclnorm << nl;

                Float ConstantPFCech = 1. / std::sqrt(eigenvalue);

                if( k == 0 ) {

                    Float Cgrad = Antihor(0,1) * maxdiameter + ConstantPFCech * Diffver(0,0) * Antihor(0,1);

                    LOG << nl;
                    LOG << "Cech PF Constant:          " << ConstantPFCech << nl;
                    LOG << "Gradient PF Constant:      " << Cgrad << nl;
                    LOG << "True Gradient PF Constant: " << 1. / Constants::pi << nl;
                    LOG << "Gradient Antihor:          " << Antihor(0,1) << nl;
                    LOG << "Gradient Diffver:          " << Diffver(0,0) << nl;
                    LOG << nl;

                }    

                if( k == 1 ) {
                    
                    Float Ccurl 
                    = Antihor(0,2) * maxdiameter
                      + comparison_quotient * std::sqrt( (n+1) * 3./2. ) * height_quotient * ConstantPFCech * Diffver(1,0) * Antihor(1,1) * Diffver(0,1) * Antihor(0,2)
                      + ( (n+1) * 3./2. ) * comparison_quotient * ( 1 + height_quotient * Antihor(1,1) ) * Diffver(0,1) * Antihor(0,2) * maxdiameter;

                    LOG << nl;
                    LOG << "Cech PF Constant:      " << ConstantPFCech << nl;
                    LOG << "Curl PF Constant:      " << Ccurl << nl;
                    LOG << "True Curl PF Constant: " << 1. / ( std::sqrt(2.) * Constants::pi ) << nl;
                    LOG << nl;

                }
            }
            
        }

        contable_results << nl;

        contable_results.lg();
        contable_parameters.lg();

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        

    } 

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
