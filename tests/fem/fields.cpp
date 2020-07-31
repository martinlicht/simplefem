

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/finitediff.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Fields and Finite Differences" << endl;
        
        {
            
            cout << "Two-dimensional setting" << endl;
        
            auto scalar = []( const FloatVector& point ) -> FloatVector { 
                return FloatVector({ 
                        2 * power_numerical(point[0],3.0) + 5 * power_numerical(point[1],2.0) + 4.
                        //std::sin( 5. * point[0] ) 
                       });
            };
            
            AlternatingForm foo( 2, 0, scalar );
            
            Float h = 10e10 * std::numeric_limits<double>::epsilon();
            
            // it seems that the difference should be taken 
            // relative to the magnitude of the coordinate
            
            AlternatingForm bar = foo.laplacian(h);
            
            FloatVector point = {1.*3.14159/2,0.};

//             std::cout << h/2.0 << space << ( foo( {0.+h,0.} ) - foo( {0.,0.} ) )[0] / h << std::endl << std::endl;
//             std::cout << h/2.0 << space << ( 2 * foo( {0.,0.} ) - foo( {0.+h,0.} ) - foo( {0.-h,0.} ) )[0] / (h*h) << std::endl;

            std::cout << foo( point ) << space << foo.exteriorderivative(h)(point) << space << bar( point ) << std::endl;

            std::cout << foo.exteriorderivative(h).exteriorderivative(h)( {2., 3.} ) << std::endl;
        
        }
        
        
        {
            
            cout << "Three-dimensional setting" << endl;
        
            auto scalarfield = []( const FloatVector& point ) -> FloatVector { 
                return FloatVector({ 
                        2 * power_numerical(point[0],3.0) + 5 * power_numerical(point[1],2.0) + point[2]
                        //std::sin( 5. * point[0] ) 
                       });
            };
            
            auto vectorfield = []( const FloatVector& point ) -> FloatVector { 
                return FloatVector({ 
                        std::sin( 1. * point[0] * point[1] * point[2] ), 
                        std::sin( 2. * point[0] * point[1] * point[2] ), 
                        std::sin( 5. * point[0] * point[1] * point[2] ), 
                       });
            };
            
            AlternatingForm form0( 3, 0, scalarfield );
            AlternatingForm form1( 3, 1, vectorfield );
            
            Float h = 10e10 * std::numeric_limits<double>::epsilon();
            
            auto d_form0  =   form0.exteriorderivative( h );
            auto dd_form0 = d_form0.exteriorderivative( h );

            auto d_form1  =   form1.exteriorderivative( h );
            auto dd_form1 = d_form1.exteriorderivative( h );

            for( int k = 0; k < 5; k++ ){
                FloatVector v(3); 
                v.random();
                cout << d_form1(v) << endl;
            }
            
            
        
        }

        {
            
            cout << "Six-dimensional setting" << endl;
        
            auto scalarfield = []( const FloatVector& point ) -> FloatVector { 
                return FloatVector({ 
                        2 * power_numerical(point[0],3.0) + 5 * power_numerical(point[1],2.0) + point[2]
                        //std::sin( 5. * point[0] ) 
                       });
            };
            
            
            AlternatingForm form0( 6, 0, scalarfield );
            
            Float h = 10e10 * std::numeric_limits<double>::epsilon();
            
            auto d_form0  =   form0.exteriorderivative( h );
            auto dd_form0 = d_form0.exteriorderivative( h );

            for( int k = 0; k < 5; k++ ){
                FloatVector v(6); 
                v.random();
                cout << dd_form0(v) << endl;
            }
            
            
        
        }

        return 0;
}
