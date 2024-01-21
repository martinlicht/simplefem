#ifndef INCLUDEGUARD_DENSE_SCALARFUNCTIONS_HPP
#define INCLUDEGUARD_DENSE_SCALARFUNCTIONS_HPP

#include "../basic.hpp"

#include "densematrix.hpp"


// matrix trace 

Float MatrixTrace( const DenseMatrix& x );

// Gerschgorin circles : row/column 

DenseMatrix Gerschgorin( const DenseMatrix& );
DenseMatrix GerschgorinRow( const DenseMatrix& );
DenseMatrix GerschgorinColumn( const DenseMatrix& );

// Matrix norms 

Float NormL1( const DenseMatrix& );
Float NormFrobenius( const DenseMatrix& );
Float NormMax( const DenseMatrix& );
Float NormLp( const DenseMatrix&, Float p );
Float NormRowCol( const DenseMatrix&, Float p, Float q );
Float NormColRow( const DenseMatrix&, Float p, Float q );

Float NormOperatorL1( const DenseMatrix& );
Float NormOperatorMax( const DenseMatrix& );

// Crude eigenvalue estimate 

Float EigenvalueEstimate( const DenseMatrix& );



#endif
