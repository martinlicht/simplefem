

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.avgsullivan.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (3D) degree elevation of interpolation preserves mass" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.automatic_dirichlet_flags();
        
        M.check();
        
        
        
        
        
        const int r_min = 1;
        
        const int r_max = 3;
        
        const int l_min = 2;
        
        const int l_max = 4;
        
        const int number_of_samples = 3;
        
           
        
        Float errors[ M.getinnerdimension()+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];

        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int k = 0;     k <= M.getinnerdimension(); k++ ) 
            for( int r = r_min; r <= r_max;                 r++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

                LOG << "Form degree: " << space << k << nl;

                LOG << "assemble mass matrices..." << nl;
                
                SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
                
                SparseMatrix inclusion  = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
                
                SparseMatrix averaging  = FEECSullivanAveragingMatrix( M, M.getinnerdimension(), k, r );
                
                errors[k][ l-l_min ][ r-r_min ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = inclusion.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );
                    
                    const auto included = inclusion * field;
                    
                    const auto averaged = averaging * included;
                    
                    const auto error_eucl = ( averaged - field ).norm();
                    
                    const auto error_mass = ( included - inclusion * averaged ).norm(massmatrix);

                    LOG << error_eucl << space << error_mass << nl; 
                    
                    // assert( error_eucl < 10e-14 and error_mass < 10e-14 ); 
                    
                    Float error = error_mass;

                    error = maximum( error, error_eucl );
                    
                    errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error );
                    
                }
                
            }
            
            if( l != l_max )
            {
                LOG << "Refinement..." << nl;
            
                M.uniformrefinement();
            }
            
        } 
    
        LOG << "Convergence tables" << nl;
    
    
        ConvergenceTable contables[ M.getinnerdimension()+1 ];
        
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
            contables[k].table_name = "Rounding errors D3K" + std::to_string(k);
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
            contables[k] << ( "R" + std::to_string(r) );

        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
                contables[k] << errors[k][l-l_min][r-r_min];
            
            contables[k] << nl; 
            
        }
        
        LOG << "Check that differences are small" << nl;
        
        for( int k      =     0; k <= M.getinnerdimension(); k++ ) 
        {

            contables[k].lg(); 
            LOG << "-------------------" << nl;

            for( int l      = l_min; l <=                 l_max; l++ ) 
            for( int r      = r_min; r <=                 r_max; r++ ) 
            {
                Assert( errors[k][l-l_min][r-r_min] < 10e-14, errors[k][l-l_min][r-r_min] );
            }

        }
        
        
        // for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        // {
        // }
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
