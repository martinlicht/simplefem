#ifndef INCLUDEGUARD_MATRIXMARKET
#define INCLUDEGUARD_MATRIXMARKET



#include <vector>
#include <limits>
#include <iostream>
#include <istream>
#include <string>

#include "../basic.hpp"

/*****
**
**  An I/O module to read matrices available in the matrix market format.
**  Only a subset of two different matrix types is supported:
**  dense and sparse matrices with real entries and general layout.
**  Other formats will not be used in our applications. 
**
******/

namespace MatrixMarket
{

  namespace Consts
  {
      
      std::string str_magicstring = "%%MatrixMarket";

      std::string str_matrix     = "matrix";
      
      std::string str_array      = "array";
      std::string str_dense      = "array";
      std::string str_coordinate = "coordinate";
      std::string str_sparse     = "coordinate";
      
      std::string str_complex    = "complex";
      std::string str_real       = "real";
      std::string str_integer    = "integer";
      std::string str_pattern    = "pattern";
      
      std::string str_general    = "general";
      std::string str_symmetric  = "symmetric";
      std::string str_skew       = "skew-symmetric";
      std::string str_hermitian  = "hermitian";
      
  }


  class Banner
  {
      
      public:
          
          Banner();

          void check() const;
          void print( std::ostream& ) const;
          
          enum class objecttype    { undefined, matrix };
          enum class matrixformat  { undefined, sparse, dense };
          enum class entrytype     { undefined, complex, real, pattern, integer };
          enum class matrixfeature { undefined, general, symmetric, skew, hermitian };
          
          bool isvalid() const;
          bool issupported() const;
          bool isdense() const;
          bool issparse() const;
          
          bool read( std::istream& );
          
      private:

          objecttype    myobjecttype;
          matrixformat  mymatrixformat;
          entrytype     myentrytype;
          matrixfeature mymatrixfeature;
          
  };



  Banner      ReadBanner( std::istream& );

  
  
//   Banner      ReadDenseMatrixBanner ( std::istream& );
//   void        ReadDenseMatrixHeader ( const std::istream&, int& numrows, int& numcols );
//   DenseMatrix ReadDenseMatrixData   ( const std::istream& );
// 
//   void WriteDenseMatrixBanner( std::ostream& );
//   void WriteDenseMatrixHeader( std::ostream&, const DenseMatrix& );
//   void WriteDenseMatrixData  ( std::ostream&, const DenseMatrix& );
// 
//   void WriteComment();
  
  
  
  
}
  
  
#endif