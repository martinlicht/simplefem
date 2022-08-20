
#include "readwrite.hpp"

#include <cassert>
#include <istream>
#include <ostream>




void writeSparseMatrix( const SparseMatrix& mat, std::ostream& output )
{
  writeSparseMatrixHeader( mat, output );
  writeSparseMatrixData  ( mat, output );
}

void writeSparseMatrixHeader( const SparseMatrix& mat, std::ostream& output )
{
  output << mat.getdimout() << tab << mat.getdimin() << nl << mat.getnumberofentries();
}

void writeSparseMatrixData  ( const SparseMatrix& mat, std::ostream& output )
{
  for( int i = 0; i < mat.getnumberofentries(); i++ )
  {
    output << mat.getentries()[i].row 
           << space 
           << mat.getentries()[i].column 
           << space 
           << mat.getentries()[i].value
           << nl;
  }
}


 

void readSparseMatrix( SparseMatrix& mat, std::istream& input )
{
  readSparseMatrixHeader( mat, input );
  readSparseMatrixData  ( mat, input );
}

void readSparseMatrixHeader( SparseMatrix& mat, std::istream& input )
{
  int r_temp;
  int c_temp;
  int n_temp;
  input >> r_temp >> c_temp >> n_temp;
  assert( r_temp == mat.getdimout() );
  assert( c_temp == mat.getdimin () );
  assert( n_temp == mat.getnumberofentries() );
}

void readSparseMatrixData  ( SparseMatrix& mat, std::istream& input )
{
  for( int i = 0; i < mat.getnumberofentries(); i++ )
  {
    SparseMatrix::MatrixEntry& e = mat.getentry(i); 
    input >> e.row >> e.column >> e.value;
  } 
}




SparseMatrix createSparseMatrixFromStream        ( std::istream& input )
{
  SparseMatrix mat = createSparseMatrixFromStreamByHeader( input );
  readSparseMatrixData( mat, input );
  return mat;
}

SparseMatrix createSparseMatrixFromStreamByHeader( std::istream& input )
{
  int r_temp;
  int c_temp;
  int n_temp;
  input >> r_temp >> c_temp >> n_temp;
  SparseMatrix mat = SparseMatrix( r_temp, c_temp, n_temp );
  return mat;
}






