

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (?D) Compare different matrix-matrix products" << nl;
        
        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial1D M1 = UnitInterval1D();
        MeshSimplicial2D M2 = UnitTriangle2D(); //StandardSquare2D_simple();
        MeshSimplicial3D M3 = UnitSimplex3D();
        
        M1.check();
        M2.check();
        M3.check();
        
        Mesh* Ms[3] = { &M1, &M2, &M3 };
        
        
        
        const int number_of_samples = 50;
        
        const int number_of_comparisons = 6;
        
        const int l_min = 0;
        
        const int l_max = 3;
        
        Float errors[l_max-l_min+1][3][number_of_comparisons];
        
        
        
        for( int l = 0; l < l_min; l++ ) {
            M1.uniformrefinement();
            M2.uniformrefinement();
            M3.uniformrefinement();            
        }
        
        
        
        for( int l = l_min; l <= l_max; l++ )
        {
        
            for( int d = 0; d < 3; d++ )
            {
                
                int m = l - l_min;
                
                for( int t = 0; t < number_of_comparisons; t++ )
                    errors[m][d][t] = -10.;
                
                Mesh& M = *(Ms[d]);
                
                
                LOG << "DIMENSION " << d+1 << " AT LEVEL " << l << nl;
        
                LOG << "...basic FEEC matrices" << nl;
        
                auto feec_vectormass = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, 0 );
                
                auto feec_diff = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, 1 );

                auto feec_diff_t = feec_diff.getTranspose();

                auto feec_inc = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, 1 );
                
                auto feec_inc_t = feec_inc.getTranspose();
            

                auto feec_vectormass_csr = MatrixCSR( feec_vectormass );
                
                auto feec_diff_csr = MatrixCSR( feec_diff );

                auto feec_diff_t_csr = MatrixCSR( feec_diff_t );

                auto feec_inc_csr = MatrixCSR( feec_inc );
                
                auto feec_inc_t_csr = MatrixCSR( feec_inc_t );
            

                LOG << "...composed FEEC matrices" << nl;
                
                auto stiffness_cpp  = feec_inc_t * feec_diff_t * feec_vectormass * feec_diff * feec_inc;
                auto stiffness_coo1 = feec_inc_t & feec_diff_t & feec_vectormass & feec_diff & feec_inc;
                auto stiffness_coo2 = ( feec_inc_t & feec_diff_t ) & feec_vectormass & ( feec_diff & feec_inc );
                auto stiffness_coo3 = feec_inc_t & ( feec_diff_t & feec_vectormass & feec_diff ) & feec_inc;

                auto stiffness_csr1 = feec_inc_t_csr & feec_diff_t_csr & feec_vectormass_csr & feec_diff_csr & feec_inc_csr;
                auto stiffness_csr2 = ( feec_inc_t_csr & feec_diff_t_csr ) & feec_vectormass_csr & ( feec_diff_csr & feec_inc_csr );
                auto stiffness_csr3 = feec_inc_t_csr & ( feec_diff_t_csr & feec_vectormass_csr & feec_diff_csr ) & feec_inc_csr;

                
                LOG << "...COMPARISONS" << nl;
                    
                for( int n = 0; n < number_of_samples; n++ ){
                    auto vec = stiffness_cpp.createinputvector();
                    vec.zero();
                    vec.random();
                    vec.normalize();
                    assert( vec.isfinite() );
                    
                    // vs COO1
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_coo1 ) * vec ).norm();
                    
                        errors[m][d][0] = maximum( vec_error, errors[m][d][0] );
                    }
                    
                    /* vs COO2 */
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_coo2 ) * vec ).norm();
                    
                        errors[m][d][1] = maximum( vec_error, errors[m][d][1] );
                    }
                    
                    /* vs COO3 */
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_coo3 ) * vec ).norm();
                    
                        errors[m][d][2] = maximum( vec_error, errors[m][d][1] );
                    }
                    
                    /* vs CSR1 */
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_csr1 ) * vec ).norm();
                    
                        errors[m][d][3] = maximum( vec_error, errors[m][d][2] );
                    }
                    
                    /* vs CSR2 */
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_csr2 ) * vec ).norm();
                    
                        errors[m][d][4] = maximum( vec_error, errors[m][d][3] );
                    }
                    
                    /* vs CSR3 */
                    {
                        auto vec_error = ( ( stiffness_cpp - stiffness_csr3 ) * vec ).norm();
                    
                        errors[m][d][5] = maximum( vec_error, errors[m][d][3] );
                    }
                    
                }
                
            } // loop over d  
        
                
            if( l != l_max ) { 

                LOG << "Refinement..." << nl;
        
                M1.uniformrefinement();
                M2.uniformrefinement();
                M3.uniformrefinement();

            }
        
        }
        
        
        
        {
            
            LOG << "Convergence tables, final results" << nl;

            ConvergenceTable contable[3];
            
            
            
            for( int d = 0; d <            3; d++ )
            {
                contable[d].table_name = "Rounding errors, D" + std::to_string(d+1);
                contable[d] << "COO1";           // 0
                contable[d] << "COO2";        // 1
                contable[d] << "COO3";        // 2 
                contable[d] << "CSR1";     // 3
                contable[d] << "CSR2";       // 4
                contable[d] << "CSR3";       // 5
                // contable[d] << "stiff comp";    // 6
                // contable[d] << "br mass";       // 7
                // contable[d] << "br stiff";      // 8
                
                
                for( int m = 0; m <= l_max-l_min; m++ ) 
                {
                    
                    for( int t = 0; t < number_of_comparisons; t++ )
                    {
                        contable[d] << errors[m][d][t];
                    }
                    
                    contable[d] << nl; 
                    
                }    
            }

            for( int d = 0; d < 3; d++ ) {
                LOG << "Dimension: " << d+1 << '\n';
                contable[d].lg();
                LOG << "----------------------------------" << nl;
            }

        }
            
            
        LOG << "Check that differences are small" << nl;
        
        for( int l = l_min; l <= l_max; l++ ) 
        for( int d = 0; d < 3; d++ )
        for( int t = 0; t < number_of_comparisons; t++ )
        {
            if( not ( errors[l-l_min][d][t] < 100 * machine_epsilon ) )
                LOG << l << space << d << space << t << space << errors[l-l_min][d][t] << nl;
            assert( errors[l-l_min][d][t] < 100 * machine_epsilon );
        }
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
