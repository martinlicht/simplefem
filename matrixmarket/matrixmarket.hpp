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
        
        const std::string str_magicstring = "%%MatrixMarket";

        const std::string str_matrix     = "matrix";
        
        const std::string str_array      = "array";
        const std::string str_dense      = "array";
        const std::string str_coordinate = "coordinate";
        const std::string str_sparse     = "coordinate";
        
        const std::string str_complex    = "complex";
        const std::string str_real       = "real";
        const std::string str_integer    = "integer";
        const std::string str_pattern    = "pattern";
        
        const std::string str_general    = "general";
        const std::string str_symmetric  = "symmetric";
        const std::string str_skew       = "skew-symmetric";
        const std::string str_hermitian  = "hermitian";
        
    }


    enum class objecttype    { undefined, matrix };
    enum class matrixformat  { undefined, sparse, dense };
    enum class entrytype     { undefined, complex, real, pattern, integer };
    enum class matrixfeature { undefined, general, symmetric, skew, hermitian };
    
    
    
    struct MatrixMarketEntry
    {
        int row;
        int column;
        double value;
    };
    
    void Read( std::istream& input, int& rows, int& columns, std::vector<MatrixMarketEntry>& entries );
    
    void WriteSparse( std::ostream& output, const int& rows, const int& columns, const std::vector<MatrixMarketEntry>& entries, const std::string& comment = "" );
    
}
  
  
#endif
