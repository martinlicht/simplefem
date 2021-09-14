
/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"



using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Simple operators" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        {
        
                LOG << "Unit Test for Scaling Component" << endl;
        
                FloatVector a(5);
                for( int i = 0; i < 5; i++ )
                        a.setentry( i, i+1 );
                
                ScalingOperator S( 5, 3.14159 );
                
                LOG << "We start with this Vector:" << endl;
                LOG << a << endl;
                
                LOG << "Scaled with PI:" << endl;
                LOG << S * a << endl;
                
                LOG << "The Scaling is " << S.getscaling() << endl;
                S.setscaling( 2.718 );
                LOG << "Now the Scaling is " << S.getscaling() << endl;
                LOG << "Accordingly:" << endl;
                LOG << S * a << endl;
                
                LOG << "Product of scaling by 11 and then 4" << endl;
                LOG << ScalingOperator(10,4) * ScalingOperator(10,11) << endl;
                LOG << ScalingOperator(10,11) * ScalingOperator(10,4) << endl;
                
                FloatVector b(10);
                for( int i = 0; i < 10; i++ )
                        b.setentry( i, i+1 );
                
                LOG << ScalingOperator(10,11) * ScalingOperator(10,4) * b << endl;
        
        }


        {
        
                LOG << "Unit Test for Diagonal Component" << endl;
        
                int dim = 10;
        
                FloatVector dia(dim);
                for( int i = 0; i < dim; i++ )
                        dia.setentry( i, i * i );
                
                DiagonalOperator D( dia );
                
                LOG << "We start with these entries:" << endl;
                LOG << dia << endl;
                
                LOG << "This is the diagonal matrix:" << endl;
                LOG << D << endl;
                
                LOG << "Pick a vector:" << endl;
                FloatVector x(dim);
                for( int i = 0; i < dim; i++ )
                        x.setentry( i, 3.01 );
                LOG << x << endl;
                LOG << "Apply the diagonal operator:" << endl;
                LOG << D * x << endl;
                
                LOG << "Now the product of diagonal operators: " << endl;
                LOG << D * D << endl;
                
                LOG << "Finished Unit Test" << endl;

                return 0;
        
        }
        
        std::clog << "Finished Unit Test: " << TestName << endl;

        return 0;
}
