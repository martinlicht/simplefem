
#include "readwrite.hpp"






void writeDenseMatrix( const DenseMatrix& mat, std::ostream& output )
{
  writeDenseMatrixHeader( mat, output );
  writeDenseMatrixData  ( mat, output );
}

void writeDenseMatrixHeader( const DenseMatrix& mat, std::ostream& output )
{
  output << mat.getdimout() << tab << mat.getdimin() << nl;
}

void writeDenseMatrixData  ( const DenseMatrix& mat, std::ostream& output )
{
  for( int r = 0; r < mat.getdimout(); r++ )
  {
    
    for( int c = 0; c < mat.getdimin (); c++ )
    {
      output << mat( r, c ) << space;
    }
    
    output << nl;
  }
}




void readDenseMatrix( DenseMatrix& mat, std::istream& input )
{
  readDenseMatrixHeader( mat, input );
  readDenseMatrixData  ( mat, input );
}

void readDenseMatrixHeader( DenseMatrix& mat, std::istream& input )
{
  int r_temp;
  int c_temp;
  input >> r_temp >> c_temp;
  assert( r_temp == mat.getdimout() );
  assert( c_temp == mat.getdimin () );
}

void readDenseMatrixData  ( DenseMatrix& mat, std::istream& input )
{
  for( int r = 0; r < mat.getdimout(); r++ )
  for( int c = 0; c < mat.getdimin (); c++ )
    input >> mat( r, c ); 
}






DenseMatrix createDenseMatrixFromStream        ( std::istream& input )
{
  DenseMatrix mat = createDenseMatrixFromStreamByHeader( input );
  readDenseMatrixData( mat, input );
  return mat;
}

DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input )
{
  int r_temp;
  int c_temp;
  input >> r_temp >> c_temp;
  DenseMatrix mat = DenseMatrix( r_temp, c_temp );
  return mat;
}

DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input, Float initval )
{
  return createDenseMatrixFromStreamByHeader( input, [initval](int,int)-> Float { return initval; } );
}

DenseMatrix createDenseMatrixFromStreamByHeader( std::istream& input, std::function<Float(int,int)> generator )
{
  DenseMatrix M = createDenseMatrixFromStreamByHeader( input );
  for( int r = 0; r < M.getdimout(); r++ )
  for( int c = 0; c < M.getdimin (); c++ )
    M( r, c ) = generator( r, c );
  return M;
}




