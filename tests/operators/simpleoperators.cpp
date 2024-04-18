

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"



int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simple operator" << nl;
    
    {
    
        LOG << "Unit Test for Scaling Operator" << nl;

        FloatVector a(5);
        for( int i = 0; i < 5; i++ )
                a.setentry( i, i+1 );
        
        ScalingOperator S( 5, Constants::pi );
        
        LOG << "We start with this Vector:" << nl;
        LOG << a << nl;
        
        LOG << "Scaled with PI:" << nl;
        const auto b = S * a;
        LOG << b << nl;
        assert( b == Constants::pi * a );
        
        LOG << "The Scaling is " << S.getscaling() << nl;
        S.setscaling( Constants::euler );
        LOG << "Now the Scaling is " << S.getscaling() << nl;
        LOG << "Accordingly:" << nl;
        const auto c = S * a;
        LOG << c << nl;
        assert( c == Constants::euler * a );
        
        LOG << "Product of scaling by 11 and then 4" << nl;
        LOG << ScalingOperator(10,4) * ScalingOperator(10,11) << nl;
        LOG << ScalingOperator(10,11) * ScalingOperator(10,4) << nl;
        
        FloatVector d(10);
        for( int i = 0; i < 10; i++ )
                d.setentry( i, i+1 );
        
        const FloatVector e = ScalingOperator(10,11) * ScalingOperator(10,4) * d;

        assert( e == 44 * d );
        
        LOG << e << nl;
    
    }


    {
    
        LOG << "Unit Test for Diagonal Operator" << nl;

        int dim = 10;

        FloatVector dia(dim);
        for( int i = 0; i < dim; i++ )
                dia.setentry( i, i * i );
        
        DiagonalOperator D( dia );
        
        LOG << "We start with these entries:" << nl;
        LOG << dia << nl;
        
        LOG << "This is the diagonal matrix:" << nl;
        LOG << D << nl;
        
        LOG << "Pick a vector:" << nl;
        FloatVector x(dim);
        for( int i = 0; i < dim; i++ )
                x.setentry( i, 3.01 );
        LOG << x << nl;

        LOG << "Apply the diagonal operator:" << nl;
        const auto y = D * x;
        LOG << D * x << nl;
        for( int i = 1; i < y.getdimension(); i++ ) assert( y[i] == i * i * x[i] );

        
        LOG << "Now the product of diagonal operators: " << nl;
        const auto F = D * D;
        LOG << F << nl;
        const auto z = F * x;
        for( int i = 1; i < z.getdimension(); i++ ) assert( z[i] == i * i * i * i * x[i] );

    }
    
    {
    
        LOG << "Unit Test for Lambda Operator" << nl;

        int dim = 10;

        FloatVector vec(dim);
        for( int i = 0; i < dim; i++ )
                vec.setentry( i, i * i );
        
        LambdaOperator L( dim,
                [](const FloatVector& input ) -> FloatVector
                {
                        return 3. * input;
                }
        );
        
        LOG << "Pick a vector:" << nl;
        FloatVector x(dim);
        for( int i = 0; i < dim; i++ )
                x.setentry( i, i );
        LOG << x << nl;
        LOG << "Apply the Lambda operator:" << nl;
        LOG << L * x << nl;

        LOG << "Apply the diagonal operator:" << nl;
        const auto y = L * x;
        LOG << L * x << nl;
        for( int i = 1; i < y.getdimension(); i++ ) assert( y[i] == 3. * x[i] );
            
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
