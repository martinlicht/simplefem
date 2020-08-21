#ifndef INCLUDEGUARD_DENSE_READWRITE
#define INCLUDEGUARD_DENSE_READWRITE

#include <iostream>
#include <functional>

#include "../basic.hpp"

#include "densematrix.hpp"

/**** 
   ** 
   ** Dense matrix input and output 
   ** - Header and Data, and a wrapper of both
   ** - own format 
   ** - 
   ** 
   *****/ 
 
void writeDenseMatrix( const DenseMatrix& mat, std::ostream& output );
void writeDenseMatrixHeader( const DenseMatrix& mat, std::ostream& output );
void writeDenseMatrixData  ( const DenseMatrix& mat, std::ostream& output );

void readDenseMatrix( DenseMatrix& mat, std::istream& input );
void readDenseMatrixHeader( DenseMatrix& mat, std::istream& input );
void readDenseMatrixData  ( DenseMatrix& mat, std::istream& input );

inline std::ostream& operator<<( std::ostream& output, const DenseMatrix& mat )
{
  writeDenseMatrix( mat, output );
  return output;
}

inline std::istream& operator>>( std::istream&  input,       DenseMatrix& mat )
{
  readDenseMatrix( mat, input );
  return input;
}




DenseMatrix createDenseMatrixFromStream        ( std::istream& input );
DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input );
DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input, Float initval );
DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input, std::function<Float(int,int)> generator );








#endif
