
#include "../../basic.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/simplesolver.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Tests for simple solvers" << nl;
    
    for( int dim = 0; dim <     10; dim++ )
    for( int t   = 0;   t <= dim+1; t++   )
    {
        DenseMatrix A(dim,dim);
        A.randommatrix();
        
        FloatVector b(dim);
        b.random();

        {
            const DenseMatrix Part = DiagonalPart(A);
            DenseMatrix PartInv = DiagonalInverse( Part );

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
            {
                if( r != c ) 
                    assert( Part(r,c) == 0. );
                else 
                    assert( Part(r,c) == A(r,c) );
            }

            assert( (Part*PartInv).is_numerically_identity() );
            assert( (PartInv*Part).is_numerically_identity() );

            auto PartInv2 = Part;
            InvertDiagonal( PartInv2 );
            Assert( ( PartInv - PartInv2 ).is_numerically_small(), A, b, Part, PartInv, PartInv2 );

            Float det = DiagonalDeterminant( Part );
            assert( is_numerically_close( det, Determinant(Part) ) );

            FloatVector x( dim );
            DiagonalSolve( Part, x, b );
            FloatVector d = Part * x - b;
            Assert( (Part*x-b).is_numerically_small(), d, Part, x, b );
        }

        {
            const DenseMatrix Part = LowerTriangularPart(A);
            DenseMatrix PartInv = LowerTriangularInverse( Part );

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
            {
                if( r < c ) 
                    assert( Part(r,c) == 0. );
                else 
                    assert( Part(r,c) == A(r,c) );
            }

            assert( (Part*PartInv).is_numerically_identity() );
            assert( (PartInv*Part).is_numerically_identity() );

            auto PartInv2 = Part;
            InvertLowerTriangular( PartInv2 );
            assert( ( PartInv - PartInv2 ).is_numerically_small() );

            Float det = LowerTriangularDeterminant( Part );
            assert( is_numerically_close( det, Determinant(Part) ) );

            FloatVector x( dim );
            LowerTriangularSolve( Part, x, b );
            FloatVector d = Part * x - b;
            Assert( (Part*x-b).is_numerically_small(), d, Part, x, b );
        }

        {
            const DenseMatrix Part = UpperTriangularPart(A);
            DenseMatrix PartInv = UpperTriangularInverse( Part );

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
            {
                if( r > c ) 
                    assert( Part(r,c) == 0. );
                else 
                    assert( Part(r,c) == A(r,c) );
            }

            Assert( (Part*PartInv).is_numerically_identity(), Part, PartInv );
            assert( (PartInv*Part).is_numerically_identity() );

            auto PartInv2 = Part;
            InvertUpperTriangular( PartInv2 );
            assert( ( PartInv - PartInv2 ).is_numerically_small() );

            Float det = UpperTriangularDeterminant( Part );
            assert( is_numerically_close( det, Determinant(Part) ) );

            FloatVector x( dim );
            UpperTriangularSolve( Part, x, b );
            FloatVector d = Part * x - b;
            Assert( (Part*x-b).is_numerically_small(), d, Part, x, b );
        }

        {
            const DenseMatrix Part = LowerUnitTriangularPart(A);
            DenseMatrix PartInv = LowerUnitTriangularInverse( Part, true );

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
            {
                if( r < c ) 
                    assert( Part(r,c) == 0. );
                else if( r == c ) 
                    assert( Part(r,c) == 1. );
                else 
                    assert( Part(r,c) == A(r,c) );
            }

            Assert( (Part*PartInv).is_numerically_identity(), Part, PartInv );
            assert( (PartInv*Part).is_numerically_identity() );

            auto PartInv2 = Part;
            InvertLowerUnitTriangular( PartInv2 );
            assert( ( PartInv - PartInv2 ).is_numerically_small() );

            FloatVector x( dim );
            LowerUnitTriangularSolve( Part, x, b );
            FloatVector d = Part * x - b;
            Assert( (Part*x-b).is_numerically_small(), d, Part, x, b );
        }

        {
            const DenseMatrix Part = UpperUnitTriangularPart(A);
            DenseMatrix PartInv = UpperUnitTriangularInverse( Part, true );

            for( int r = 0; r < dim; r++ )
            for( int c = 0; c < dim; c++ )
            {
                if( r > c ) 
                    assert( Part(r,c) == 0. );
                else if( r == c ) 
                    assert( Part(r,c) == 1. );
                else 
                    assert( Part(r,c) == A(r,c) );
            }

            assert( (Part*PartInv).is_numerically_identity() );
            assert( (PartInv*Part).is_numerically_identity() );

            auto PartInv2 = Part;
            InvertUpperUnitTriangular( PartInv2 );
            assert( ( PartInv - PartInv2 ).is_numerically_small() );

            FloatVector x( dim );
            UpperUnitTriangularSolve( Part, x, b );
            FloatVector d = Part * x - b;
            Assert( (Part*x-b).is_numerically_small(), d, Part, x, b );
        }
        
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
