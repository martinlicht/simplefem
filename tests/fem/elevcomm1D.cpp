

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (1D) degree elevations commute" << endl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << endl;
        
        auto M = UnitInterval1D();
        
        M.check();
        


        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 4;
        
        const int number_of_samples = 3;
        
           
        
        Float errors[ M.getinnerdimension()+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << endl;
            
            for( int k = 0; k <= M.getinnerdimension(); k++ ) 
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << endl;

                LOG << "Form degree: " << space << k << endl;

                LOG << "assemble matrices..." << endl;
        
                SparseMatrix elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 1 );
                SparseMatrix elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+1, 1 );
                SparseMatrix elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+2, 1 );
                SparseMatrix elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 3 );
                
                errors[k][ l ][ r ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = elevation_r_g.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );
                    
                    const auto path_direct   = elevation_r_g * field;
                    
                    const auto path_indirect = elevation_r_3 * elevation_r_2 * elevation_r_1 * field;
                    
                    const auto error_mass = ( path_direct - path_indirect ).norm();
                    
                    errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error_mass );
                    
                }
                
            }

            if( l != l_max )
            {
                LOG << "Refinement..." << endl;
            
                M.uniformrefinement();
            }
            
        } 
        
        
        
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable[ M.getinnerdimension()+1 ];
        
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
            contable->table_name = "Rounding errors D1K" + std::to_string(k);
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
            contable[k] << ( "R" + std::to_string(r-r_min) );

        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
                contable[k] << errors[k][l-l_min][r-r_min];
            
            contable[k] << nl; 
            
        }
        
        
        
        for( int k = 0; k <= M.getinnerdimension(); k++ ) 
        {
            contable[k].lg(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        LOG << "Check that differences are small" << nl;
        
        for( int l      = l_min; l <=                 l_max; l++ ) 
        for( int r      = r_min; r <=                 r_max; r++ ) 
        for( int k      =     0; k <= M.getinnerdimension(); k++ ) 
        {
            assert( errors[k][l-l_min][r-r_min] < 10e-14 );
        }
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
