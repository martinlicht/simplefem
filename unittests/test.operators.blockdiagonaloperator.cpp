

/**/

#include <iostream>
#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/scalingoperator.hpp"
#include "../operators/blockdiagonaloperator.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Block Diagonal Operator" << endl;
    
    FloatVector a(20,1.);
    
    ScalingOperator S1(  5, 3.14159 );
    ScalingOperator S2( 11, 2.718 );
    ScalingOperator S3(  4, 2.5029 );
    
    std::vector<LinearOperator*> ops;
    ops.push_back( &S1 );
    ops.push_back( &S2 );
    ops.push_back( &S3 );
    BlockDiagonalOperator BDO( 20, 20, ops );
    
    cout << "The following operators:" << endl;
    cout << S1 << S2 << S3 << endl;
    
    cout << "come into a block diagonal operator" << endl;
    cout << BDO << endl;
    
    cout << "Application:" << endl;
    cout << BDO * a << endl;
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
