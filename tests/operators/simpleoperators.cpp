
/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"



using namespace std;

int main()
{
        cout << "Unit Test for Simple operator" << endl;
        
        {
        
                cout << "Unit Test for Scaling Component" << endl;
        
                FloatVector a(5);
                for( int i = 0; i < 5; i++ )
                        a.setentry( i, i+1 );
                
                ScalingOperator S( 5, 3.14159 );
                
                cout << "We start with this Vector:" << endl;
                cout << a << endl;
                
                cout << "Scaled with PI:" << endl;
                cout << S * a << endl;
                
                cout << "The Scaling is " << S.getscaling() << endl;
                S.setscaling( 2.718 );
                cout << "Now the Scaling is " << S.getscaling() << endl;
                cout << "Accordingly:" << endl;
                cout << S * a << endl;
                
                cout << "Product of scaling by 11 and then 4" << endl;
                cout << ScalingOperator(10,4) * ScalingOperator(10,11) << endl;
                cout << ScalingOperator(10,11) * ScalingOperator(10,4) << endl;
                
                FloatVector b(10);
                for( int i = 0; i < 10; i++ )
                        b.setentry( i, i+1 );
                
                cout << ScalingOperator(10,11) * ScalingOperator(10,4) * b << endl;
        
        }


        {
        
                cout << "Unit Test for Diagonal Component" << endl;
        
                int dim = 10;
        
                FloatVector dia(dim);
                for( int i = 0; i < dim; i++ )
                        dia.setentry( i, i * i );
                
                DiagonalOperator D( dia );
                
                cout << "We start with these entries:" << endl;
                cout << dia << endl;
                
                cout << "This is the diagonal matrix:" << endl;
                cout << D << endl;
                
                cout << "Pick a vector:" << endl;
                FloatVector x(dim);
                for( int i = 0; i < dim; i++ )
                        x.setentry( i, 3.01 );
                cout << x << endl;
                cout << "Apply the diagonal operator:" << endl;
                cout << D * x << endl;
                
                cout << "Now the product of diagonal operators: " << endl;
                cout << D * D << endl;
                
                cout << "Finished Unit Test" << endl;

                return 0;
        
        }
        
        cout << "Finished Unit Test" << endl;

        return 0;
}
