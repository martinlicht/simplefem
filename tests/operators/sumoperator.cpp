

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/scalingoperator.hpp"
#include "../../operators/sumoperator.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for SumOperator Class" << endl;
        
        if( true ) {
          
            cout << "Test with two scaling matrices" << endl;
        
            ScalingOperator S1( 5, 3.14159 );
            
            ScalingOperator S2( 5, 4.6692 );
            
            FloatVector x( 5 );
            for( int i = 0; i < 5; i++ )
                x[i] = i*i;
        
            cout << x << endl;
            
            cout << SumOperator( S1, S2 ) * x << endl;
        
        }
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
