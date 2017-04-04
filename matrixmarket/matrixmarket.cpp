

#include "matrixmarket.hpp"

#include <iostream>
#include <algorithm>


MatrixMarket::Banner::Banner() : 
  myobjecttype(    MatrixMarket::Banner::objecttype::undefined ), 
  mymatrixformat(  MatrixMarket::Banner::matrixformat::undefined ), 
  myentrytype(     MatrixMarket::Banner::entrytype::undefined ), 
  mymatrixfeature( MatrixMarket::Banner::matrixfeature::undefined ) 
{
}

void MatrixMarket::Banner::check() const
{
}

bool MatrixMarket::Banner::isvalid() const
{
  if( myobjecttype   != MatrixMarket::Banner::objecttype::matrix )
    return false;
  
  if( mymatrixformat == MatrixMarket::Banner::matrixformat::dense 
      &&
      myentrytype    == MatrixMarket::Banner::entrytype::pattern
    )
    return false;
  
  if( mymatrixfeature == MatrixMarket::Banner::matrixfeature::hermitian
      &&
      myentrytype     == MatrixMarket::Banner::entrytype::real
    )
    return false;
  
  if( myentrytype     == MatrixMarket::Banner::entrytype::pattern
      &&
      (
        mymatrixfeature == MatrixMarket::Banner::matrixfeature::hermitian
        ||
        mymatrixfeature == MatrixMarket::Banner::matrixfeature::skew
      )
    )
    return false;
  
  return true;
}

bool MatrixMarket::Banner::issupported() const
{
  bool valid = isvalid();
  return valid
         &&
         myentrytype     == MatrixMarket::Banner::entrytype::real
         &&
         mymatrixfeature == MatrixMarket::Banner::matrixfeature::general;
         ;
}

bool MatrixMarket::Banner::isdense() const
{
  return mymatrixformat == MatrixMarket::Banner::matrixformat::dense;
}

bool MatrixMarket::Banner::issparse() const
{
  return mymatrixformat == MatrixMarket::Banner::matrixformat::sparse;
}


bool MatrixMarket::Banner::read( std::istream& input ) 
{
  
  std::string str_magicstring, str_objecttype, str_matrixformat, str_entrytype, str_matrixfeature;
  
  input >> str_magicstring;
  input >> str_objecttype;
  input >> str_matrixformat;
  input >> str_entrytype;
  input >> str_matrixfeature;
  
  std::cout << "READ: " << std::endl;
  std::cout << str_magicstring        << space
            << str_objecttype    << space
            << str_matrixformat  << space
            << str_entrytype     << space
            << str_matrixfeature << space
            << std::endl;
  
  std::transform( str_magicstring.begin(),        str_magicstring.end(),        str_magicstring.begin(),        ::tolower );
  std::transform( str_objecttype.begin(),    str_objecttype.end(),    str_objecttype.begin(),    ::tolower );
  std::transform( str_matrixformat.begin(),  str_matrixformat.end(),  str_matrixformat.begin(),  ::tolower );
  std::transform( str_entrytype.begin(),     str_entrytype.end(),     str_entrytype.begin(),     ::tolower );
  std::transform( str_matrixfeature.begin(), str_matrixfeature.end(), str_matrixfeature.begin(), ::tolower );
  
  std::cout << "TRANFORMED: " << std::endl;
  std::cout << str_magicstring        << space
            << str_objecttype    << space
            << str_matrixformat  << space
            << str_entrytype     << space
            << str_matrixfeature << space
            << std::endl;
  
  /* check: magic string is given */
  
  if( ! str_magicstring.compare( MatrixMarket::Consts::str_magicstring ) )
    return false;
  
  /* check: object type is matrix */
  if( ! str_objecttype.compare( MatrixMarket::Consts::str_matrix ) )
    return false;
  
  
  /* check: matrix format */
  if     ( ! str_matrixformat.compare( MatrixMarket::Consts::str_dense ) )
    mymatrixformat = MatrixMarket::Banner::matrixformat::dense;
  else if( ! str_matrixformat.compare( MatrixMarket::Consts::str_sparse ) )
    mymatrixformat = MatrixMarket::Banner::matrixformat::sparse;
  else
    return false;
  
  /* check: entrytype */
  if     ( ! str_entrytype.compare( MatrixMarket::Consts::str_real ) )
    myentrytype = MatrixMarket::Banner::entrytype::real;
  else if( ! str_entrytype.compare( MatrixMarket::Consts::str_complex ) )
    myentrytype = MatrixMarket::Banner::entrytype::complex;
  else if( ! str_entrytype.compare( MatrixMarket::Consts::str_integer ) )
    myentrytype = MatrixMarket::Banner::entrytype::integer;
  else if( ! str_entrytype.compare( MatrixMarket::Consts::str_pattern ) )
    myentrytype = MatrixMarket::Banner::entrytype::pattern;
  else
    return false;
  
  
  /* check: storage scheme */
  if     ( ! str_matrixfeature.compare( MatrixMarket::Consts::str_general ) )
    mymatrixfeature = MatrixMarket::Banner::matrixfeature::general;
  else if( ! str_matrixfeature.compare( MatrixMarket::Consts::str_symmetric ) )
    mymatrixfeature = MatrixMarket::Banner::matrixfeature::symmetric;
  else if( ! str_matrixfeature.compare( MatrixMarket::Consts::str_hermitian ) )
    mymatrixfeature = MatrixMarket::Banner::matrixfeature::hermitian;
  else if( ! str_matrixfeature.compare( MatrixMarket::Consts::str_skew ) )
    mymatrixfeature = MatrixMarket::Banner::matrixfeature::skew;
  else
    return false;
  
    
  return true;
}









