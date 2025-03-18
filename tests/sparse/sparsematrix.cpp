
#include "../../base/include.hpp"
#include "../../utility/random.hpp"
#include "../../operators/floatvector.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: SparseMatrix" << nl;

    {
        LOG << "1. Action of Identity Matrix" << nl;
        const int dim = 5;
        
        SparseMatrix M( dim, dim );

        for( int i = 0; i < dim; i++ ) M.appendentry( i, i, 1. );
        
        FloatVector vec( dim, 1.23 );
        for( int i = 0; i < dim; i++ ) vec[i] = i * 1.2345;

        auto Mvec = M * vec;

        for( int i = 0; i < dim; i++ ) assert( vec[i] == Mvec[i] );
        
    }
    
    {
        LOG << "2. Some sparse matrix" << nl;
        
        SparseMatrix M( 2, 3 );

        for( int i = 0; i < 5; i++ )
            for( int j = 0; j < 7; j++ )
                M.appendentry( (3*i) % 2, (2*j) % 3, i / 3. + j*j );

        LOG << "This is the content of some matrix:" << nl;
        LOG << M << nl;

        LOG << "Sort the Entries" << nl;
        M.sortentries();
        LOG << M << nl;

        M.clearentries();
        LOG << "Empty Matrix again" << nl;
        LOG << M << nl;

        LOG << "Next Matrix:" << nl;
        M.appendentry( 0, 0, 1. );
        M.appendentry( 0, 1, 2. );
        M.appendentry( 0, 2, 3. );
        M.appendentry( 1, 0, 4. );
        M.appendentry( 1, 1, 5. );
        M.appendentry( 1, 2, 6. );
        LOG << M << nl;

        FloatVector vec(3);
        vec.setentry(0,13);
        vec.setentry(1,17);
        vec.setentry(2,19);
        LOG << "Some vector:" << nl;
        LOG << vec << nl;

        LOG << "Matrix-Vector Product:" << nl;
        LOG << M * vec << nl;
    
    }
    
    if(false)
    {
        
        LOG << "3. Bulk testing" << nl;

        const int max_dimout = 20;
        const int max_dimin  = 20;
        
        for( int dimout = 1; dimout <= max_dimout; dimout++ )
        for( int dimin  = 1; dimin  <= max_dimin;  dimin ++ )
        for( int t = 0; t < 7; t++ )
        {
            FloatVector input( dimin );
            input.random();
            
            FloatVector output( dimout, 0. );

            SparseMatrix M( dimout, dimin );

            for( int s = 0; s < 3 * dimin * dimout; s++ )
            {
                int   r = random_integer() % dimout;
                int   c = random_integer() % dimin;
                Float v = 10 * random_uniform() - 5.;

                output[r] += v * input[c];

                M.appendentry( r, c, v );
            }

            auto new_output = M * input;

            assert( (output-new_output).is_numerically_small() );

            M.sortandcompressentries();

            auto compressed_output = M * input;

            assert( (output-compressed_output).is_numerically_small() );
        }

    }
    
    if(false)
    {
        
        LOG << "4. Compare sparse and CSR matrix" << nl;
        
        for( int dimout = 1; dimout <= 10; dimout++ )
        for( int dimin  = 1; dimin  <= 10; dimin ++ )
        for( int t = 0; t < 7; t++ )
        {
            FloatVector input( dimin );
            input.random();
            
            FloatVector output( dimout, 0. );

            SparseMatrix M( dimout, dimin );

            for( int s = 0; s < 3 * dimin * dimout; s++ )
            {
                int   r = random_integer() % dimout;
                int   c = random_integer() % dimin;
                Float v = 10 * random_uniform() - 5.;

                output[r] += v * input[c];

                M.appendentry( r, c, v );
            }

            auto M_output = M * input;

            assert( (output-M_output).is_numerically_small() );

            MatrixCSR A = MatrixCSR( M );

            auto A_output = A * input;

            assert( (output-A_output).is_numerically_small() );
        }

    }
    
    {
        
        LOG << "5. Output sparse and CSR matrix" << nl;
        
        SparseMatrix M( 20, 20 );

        for( int s = 0; s < 100; s++ )
        {
            int   r = random_integer() % M.getdimout();
            int   c = random_integer() % M.getdimin();
            Float v = 10 * random_uniform() - 5.;

            M.appendentry( r, c, v );
        }

        M.save_graphics( "sparsematrix.bmp" );

        MatrixCSR csr = MatrixCSR( M );

        csr.save_graphics( "csrmatrix.bmp" );
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
