

/**/

#include "../operators/diagonaloperator.hpp"

#include <iostream>
#include "../basic.hpp"
#include "../operators/floatvector.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for ScalingOperator" << endl;
        
        int dim = 10;
        
        FloatVector dia(dim);
        for( int i = 0; i < dim; i++ )
                dia.setentry( i, i * i );
        
        DiagonalOperator D( dim, dia );
        
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
