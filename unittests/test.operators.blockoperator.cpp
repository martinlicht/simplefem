

/**/

#include <iostream>
#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/scalingoperator.hpp"
#include "../operators/blockoperator.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Block Operator" << endl;
    
    FloatVector a(20,1.);
    
    ScalingOperator S1(  5,    4 );
    ScalingOperator S2(  5,    7 );
    ScalingOperator S3(  5,    5 );
    ScalingOperator S4(  5, 1000 );
    
    std::vector<std::vector<LinearOperator*>> ops;
    ops.resize(3);
    ops[0].resize(4);
    ops[0][0] = &S1;
    ops[0][1] = &S1;
    ops[0][2] = &S1;
    ops[0][3] = &S3;
    ops[1].resize(4);
    ops[1][0] = &S2;
    ops[1][1] = &S2;
    ops[1][2] = &S2;
    ops[1][3] = &S4;
    ops[2].resize(4);
    ops[2][0] = &S3;
    ops[2][1] = &S1;
    ops[2][2] = &S2;
    ops[2][3] = &S3;
    
    BlockOperator BO( 15, 20, ops );
    
    cout << "The following operators:" << endl;
    cout << S1 << S2 << S3 << S4 << endl;
    
    cout << "come into a block operator" << endl;
    cout << BO << endl;
    
    cout << "Application:" << endl;
    cout << BO * a << endl;
    
    cout << "Finished Unit Test" << endl;

    return 0;
}
