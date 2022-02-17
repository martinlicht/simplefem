#ifndef INCLUDEGUARD_SPARSE_READWRITE_HPP
#define INCLUDEGUARD_SPARSE_READWRITE_HPP

#include <istream>
#include <ostream>

#include "../basic.hpp"

#include "sparsematrix.hpp"

/**** 
   ** 
   ** Sparse matrix input and output 
   ** - Header and Data, and a wrapper of both
   ** - own format 
   ** - 
   ** 
   *****/ 
 
void writeSparseMatrix( const SparseMatrix& mat, std::ostream& output );
void writeSparseMatrixHeader( const SparseMatrix& mat, std::ostream& output );
void writeSparseMatrixData  ( const SparseMatrix& mat, std::ostream& output );

void readSparseMatrix( SparseMatrix& mat, std::istream& input );
void readSparseMatrixHeader( SparseMatrix& mat, std::istream& input );
void readSparseMatrixData  ( SparseMatrix& mat, std::istream& input );

// inline std::ostream& operator<<( std::ostream& output, const SparseMatrix& mat )
// {
//   writeSparseMatrix( mat, output );
//   return output;
// }

// inline std::istream& operator>>( std::istream&  input,       SparseMatrix& mat )
// {
//   readSparseMatrix( mat, input );
//   return input;
// }




SparseMatrix createSparseMatrixFromStream        ( std::istream& input );
SparseMatrix createSparseMatrixFromStreamByHeader( std::istream& input );









#endif
