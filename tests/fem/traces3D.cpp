

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.flags.hpp"
#include "../../utility/convergencetable.hpp"

#include "../../fem/global.trace.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
        
    LOG << "Unit Test: (3D) degree elevation of interpolation preserves mass" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitCube3D();
    
    // M.automatic_dirichlet_flags();
    
    M.check();
    
    
    const int r_min = 1;
    
    const int r_max = 3;
    
    const int l_min = 0;
    
    const int l_max = 3;
    
    const int number_of_samples = 3;
        
    
    Float errors[ M.getinnerdimension()+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];

    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        for( int k = 0;     k <= M.getinnerdimension(); k++ ) 
        for( int r = r_min; r <= r_max;                 r++ ) 
        {
            
            LOG << "Level:       " << l_min << " <= " << l << " <= " << l_max << nl;
            
            LOG << "Polydegree:  " << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "Form degree: " << k << nl;

            LOG << "Level: " << l << "/" << l_max << nl;
            LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;


            LOG << "assemble matrices..." << nl;
            
            const auto inclusion  = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            const auto trace      = FEECBrokenTraceMatrix( M, M.getinnerdimension(), k, r, false );

            const auto massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension()-1, k, r );
            
            
            
            
            errors[k][ l-l_min ][ r-r_min ] = 0.;
            
            for( int i = 0; i < number_of_samples; i++ ){

                auto field = inclusion.createinputvector();
                field.random();

                auto included_field = inclusion * field;

                auto traces_of_field = trace * included_field;

                assert( field.isfinite() );
                
                const auto error_eucl = field.norm();
                
                LOG << massmatrix.getdimout() << space << massmatrix.getdimin() << space << traces_of_field.getdimension() << nl;

                const auto error_mass = traces_of_field.norm(massmatrix);

                LOG << error_eucl << space << error_mass << nl; 
                
                Float error = error_mass;

                error = maximum( error, error_eucl );
                
                errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error );
                
            }
            
        }
        
        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();

            LOG << "Distortion..." << nl;
        
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
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
