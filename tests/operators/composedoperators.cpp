

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/scalingoperator.hpp"
#include "../../operators/composedoperators.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Produkt Operator Class" << endl;
        
        if( true ) {
          
            cout << "Test with two scaling matrices" << endl;
        
            ScalingOperator S1( 5, 3.14159 );
            ScalingOperator S2( 5, 4.6692 );
            ScalingOperator S3( 5, Constants::feigenbaum_first );
            ScalingOperator S4( 5, Constants::feigenbaum_second );
            ScalingOperator S5( 5, Constants::sirpinski );
            
            FloatVector x( 5 );
            for( int i = 0; i < 5; i++ )
                x[i] = i*i;
                
            cout << x << endl;
            
            auto test1 = S5 * ( ProduktOperator( S1, S2 ) * ProduktOperator( S3, S4 ) );
            auto test2 = S1 * S2;
            
            cout << test1 * x << endl;
            cout << test2 * x << endl;
            
        }
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
