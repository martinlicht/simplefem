

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test: (?D) Lagrange matrices agree with FEEC analogues" << endl;
        
        cout << std::setprecision(10);

        cout << "Initial mesh..." << endl;
        
        MeshSimplicial1D M1 = UnitInterval1D();
        MeshSimplicial2D M2 = UnitTriangle2D(); //StandardSquare2D_simple();
        MeshSimplicial3D M3 = UnitSimplex3D();
        
        M1.check();
        M2.check();
        M3.check();
        
        Mesh* Ms[3] = { &M1, &M2, &M3 };
        
        
        
        const int number_of_samples = 50;
        
        const int number_of_comparisons = 7;
        
        const int l_min = 0;
        
        const int l_max = 4;
        
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
                    errors[l][d][t] = 0.;
                
                Mesh& M = *(Ms[d]);
                
                
                cout << "DIMENSION " << d+1 << " AT LEVEL " << l << endl;
        
                cout << "...basic FEEC matrices" << endl;
        
                auto feec_broken_mass = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, 1 );

                auto feec_vectormass = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, 0 );
                
                auto feec_diff = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, 1 );

                auto feec_diff_t = feec_diff.getTranspose();

                auto feec_inc = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, 1 );
                
                auto feec_inc_t = feec_inc.getTranspose();
            
                assert( feec_broken_mass.isfinite() );
                assert( feec_vectormass.isfinite() );
                assert( feec_diff.isfinite() );
                assert( feec_diff_t.isfinite() );
                assert( feec_inc.isfinite() );
                assert( feec_inc_t.isfinite() );
                    
                cout << "...composed FEEC matrices" << endl;
                
                auto feec_broken_stiffness = feec_diff_t & feec_vectormass       & feec_diff;
                auto feec_stiffness        = feec_inc_t  & feec_broken_stiffness & feec_inc;
                auto feec_mass             = feec_inc_t  & feec_broken_mass      & feec_inc;
                
                assert( feec_stiffness.isfinite()        );
                assert( feec_mass.isfinite()             );
                assert( feec_broken_stiffness.isfinite() );
                
                cout << "...basic Lagrange matrices" << endl;
                
                auto lagr_broken_mass      = LagrangeBrokenMassMatrix( M, 1 );
                auto lagr_broken_stiffness = LagrangeBrokenStiffnessMatrix( M, 1 );
                
                auto lagr_mass      = LagrangeMassMatrix( M, 1 );
                auto lagr_stiffness = LagrangeStiffnessMatrix( M, 1 );
                
                auto lagr_inc   = LagrangeInclusionMatrix( M, M.getinnerdimension(), 1 );
                auto lagr_inc_t = lagr_inc.getTranspose();
                
                cout << "...composed Lagrange matrices" << endl;
                
                auto lagr_composed_mass      = lagr_inc_t & lagr_broken_mass & lagr_inc;
                auto lagr_composed_stiffness = lagr_inc_t & lagr_broken_stiffness & lagr_inc;
                
                cout << "...COMPARISONS" << nl;
                    
                for( int n = 0; n < number_of_samples; n++ ){
                    auto vec = lagr_inc.createinputvector();
                    vec.zero();
                    vec.random();
                    vec.normalize();
                    vec.isfinite();
                    
                    // inclusion matrices
                    {
                        auto vec_error = ( ( lagr_inc - feec_inc ) * vec ).norm();
                    
                        errors[m][d][0] = maximum( vec_error, errors[m][d][0] );
                    }
                    
                    /*mass matrices*/
                    {
                        auto vec_error = ( ( lagr_mass - feec_mass ) * vec ).norm();
                    
                        errors[m][d][1] = maximum( vec_error, errors[m][d][1] );
                    }
                    
                    /*mass matrices*/
                    {
                        auto vec_error = ( ( lagr_mass - lagr_composed_mass ) * vec ).norm();
                    
                        errors[m][d][2] = maximum( vec_error, errors[m][d][2] );
                    }
                    
                    /*stiffness matrices*/
                    {
                        auto vec_error = ( ( lagr_stiffness - lagr_composed_stiffness ) * vec ).norm();
                    
                        errors[m][d][3] = maximum( vec_error, errors[m][d][3] );
                    }
                    
                    /*stiffness matrices*/
                    {
                        auto vec_error = ( ( lagr_composed_stiffness - feec_stiffness ) * vec ).norm();
                    
                        errors[m][d][4] = maximum( vec_error, errors[m][d][4] );
                    }
                    
                }
                
                for( int n = 0; n < number_of_samples; n++ ){
                    auto vec = feec_broken_mass.createinputvector();
                    vec.zero();
                    vec.random();
                    vec.normalize();
                    vec.isfinite();
                    
                    /*broken mass*/
                    {
                        auto vec_error = ( ( lagr_broken_mass - feec_broken_mass ) * vec ).norm();
                    
                        errors[m][d][5] = maximum( vec_error, errors[m][d][5] );
                    }
                    
                    /*broken stiffness*/
                    {
                        auto vec_error = ( ( lagr_broken_stiffness - feec_broken_stiffness ) * vec ).norm();
                    
                        errors[m][d][6] = maximum( vec_error, errors[m][d][6] );
                    }
                    
                }
                
                
                
//                 if( d==2 and l==0 ){
//                     lagr_inc.print( std::cout );
//                     lagr_stiffness.print( std::cout );
//                     lagr_broken_stiffness.print( std::cout );
//                     feec_stiffness.print( std::cout );
//                     exit(0);
//                 }
                
                
                
            } // loop over d  
    
                
                
            cout << "Convergence tables" << nl;

            ConvergenceTable contable[3];
            
            for( int d = 0; d <        3; d++ )
            for( int m = 0; m <= l-l_min; m++ ) 
            {
                
                for( int t = 0; t < number_of_comparisons; t++ )
                {
                    contable[d] << errors[m][d][t];
                }
                
                contable[d] << nl; 
                
            }
                
            for( int d = 0; d < 3; d++ ) {
                contable[d].print( cout );
                std::cout << "----------------------------------" << std::endl;
            }
                
                
                
            cout << "Refinement..." << endl;
        
            M1.uniformrefinement();
            M2.uniformrefinement();
            M3.uniformrefinement();
        
        }
        
        
        
        {
            
            cout << "Convergence tables, final results" << nl;

            ConvergenceTable contable[3];
            
            for( int d = 0; d <            3; d++ )
            for( int m = 0; m <= l_max-l_min; m++ ) 
            {
                
                for( int t = 0; t < number_of_comparisons; t++ )
                {
                    contable[d] << errors[m][d][t];
                }
                
                contable[d] << nl; 
                
            }
                
            for( int d = 0; d < 3; d++ ) {
                contable[d].print( cout );
                std::cout << "----------------------------------" << std::endl;
            }
            
        }
            
            
        cout << "Check that differences are small" << nl;
        
        for( int l = l_min; l <= l_max; l++ ) 
        for( int d = 0; d < 3; d++ )
        for( int t = 0; t < number_of_comparisons; t++ )
        {
            assert( errors[l][d][t] < 10e-10 );
        }
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
