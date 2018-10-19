

#include "matrixmarket.hpp"

#include <iostream>
#include <algorithm>



void MatrixMarket::Read( std::istream& input, int& rows, int& columns, std::vector<MatrixMarket::MatrixMarketEntry>& entries )
{
    
    std::string str_magicstring;
    
    std::string str_objecttype;
    std::string str_matrixformat;
    std::string str_entrytype;
    std::string str_matrixfeature;
    
    objecttype    my_objecttype    = objecttype::undefined;
    matrixformat  my_matrixformat  = matrixformat::undefined;
    entrytype     my_entrytype     = entrytype::undefined;
    matrixfeature my_matrixfeature = matrixfeature::undefined;
    
    std::string my_comment = "";
  
    
    input >> str_magicstring;
    
    std::clog << str_magicstring << nl;
  
//     std::transform( str_magicstring.begin(),   str_magicstring.end(),   str_magicstring.begin(),   ::tolower );
  
    /* check: magic string is given */
    assert( 0 == str_magicstring.compare( Consts::str_magicstring ) );
    
    
    
    input >> str_objecttype;
    
    std::transform( str_objecttype.begin(),    str_objecttype.end(),    str_objecttype.begin(),    ::tolower );
    
    my_objecttype = objecttype::matrix;
    
    /* check: object type is matrix */
    
    assert( 0 == str_objecttype.compare( Consts::str_matrix ) );
  
    
    
    input >> str_matrixformat;
    
    std::transform( str_matrixformat.begin(),  str_matrixformat.end(),  str_matrixformat.begin(),  ::tolower );
    
    if( 0 == str_matrixformat.compare( Consts::str_dense  ) ) my_matrixformat = matrixformat::dense;
    if( 0 == str_matrixformat.compare( Consts::str_sparse ) ) my_matrixformat = matrixformat::sparse;
    
    /* check: matrix format is dense or sparse */
    assert( 0 == str_matrixformat.compare( Consts::str_dense ) || 0 == str_matrixformat.compare( Consts::str_sparse ) );
    
    
    
    
    
    input >> str_entrytype;
    
    std::transform( str_entrytype.begin(),     str_entrytype.end(),     str_entrytype.begin(),     ::tolower );
    
    if( 0 == str_entrytype.compare( Consts::str_real    ) ) my_entrytype = entrytype::real;
    if( 0 == str_entrytype.compare( Consts::str_complex ) ) my_entrytype = entrytype::complex;
    if( 0 == str_entrytype.compare( Consts::str_integer ) ) my_entrytype = entrytype::integer;
    if( 0 == str_entrytype.compare( Consts::str_pattern ) ) my_entrytype = entrytype::pattern;

    /* check: entrytype is real, complex, integer, or pattern */
    assert( 0 == str_entrytype.compare( Consts::str_real ) || 0 == str_entrytype.compare( Consts::str_complex ) || 0 == str_entrytype.compare( Consts::str_integer ) || 0 == str_entrytype.compare( Consts::str_pattern ) );
    
    
    
    input >> str_matrixfeature;
    
    std::transform( str_matrixfeature.begin(), str_matrixfeature.end(), str_matrixfeature.begin(), ::tolower );
  
    if( 0 == str_matrixfeature.compare( Consts::str_general   ) ) my_matrixfeature = matrixfeature::general;
    if( 0 == str_matrixfeature.compare( Consts::str_symmetric ) ) my_matrixfeature = matrixfeature::symmetric;
    if( 0 == str_matrixfeature.compare( Consts::str_hermitian ) ) my_matrixfeature = matrixfeature::hermitian;
    if( 0 == str_matrixfeature.compare( Consts::str_skew      ) ) my_matrixfeature = matrixfeature::skew;
    
    /* check: storage scheme is general, symmetric, hermitian, or skew-symmetric */
    assert( 0 == str_matrixfeature.compare( Consts::str_general ) || 0 == str_matrixfeature.compare( Consts::str_symmetric ) || 0 == str_matrixfeature.compare( Consts::str_hermitian ) || 0 == str_matrixfeature.compare( Consts::str_skew ) );
    
    
    /* check that all data are defined */
    assert( my_objecttype    != objecttype::undefined );
    assert( my_matrixformat  != matrixformat::undefined );
    assert( my_entrytype     != entrytype::undefined );
    assert( my_matrixfeature != matrixfeature::undefined );
    
    
  
    /* read the comments for what it's worth */
    
    while( input.peek() == '%' ) {
    
        input.get();
        std::string line;
        std::getline( input, line );
        my_comment = my_comment + line;
        
    }
  
  
    std::clog << "READ:" << std::endl;
    std::clog << str_magicstring   << space
                << str_objecttype    << space
                << str_matrixformat  << space
                << str_entrytype     << space
                << str_matrixfeature << space
                << std::endl;
    
    
    std::clog << "TRANFORMED:" << std::endl;
    std::clog << str_magicstring   << space
                << str_objecttype    << space
                << str_matrixformat  << space
                << str_entrytype     << space
                << str_matrixfeature << space
                << std::endl;
                
    std::clog << "COMMENTS: " << std::endl << my_comment << std::endl;
  
    
    
    /* check that only combinations of parameters have been read */
    
    assert( my_objecttype == objecttype::matrix ); // object must matrix
    assert( my_matrixformat != matrixformat::dense || my_entrytype != entrytype::pattern ); // dense matrices are not given by patterns 
    assert( my_matrixfeature != matrixfeature::hermitian || my_entrytype != entrytype::real ); // hermitian matrices must have complex entries 
    assert( my_entrytype != entrytype::pattern || !( my_matrixfeature == matrixfeature::hermitian || my_matrixfeature == matrixfeature::skew ) ); // patterns are only for general or symmetric matrices
    
    /* check that only supported data formats are used */ 
    
    assert( my_entrytype == entrytype::real || my_entrytype == entrytype::integer ); // complex and pattern matrices are not yet supported 
    assert( my_matrixfeature != matrixfeature::hermitian ); // only general matrices are currently supported
    
    // TODO: add symmetric and skew-symmetric
    
    
    
    /**
     At this point, we have read the header of the MM file.
     Next we extract the matrix dimension and the actual matrix entries. 
     ***/
    
    assert( entries.size() == 0 );
    
    if( my_matrixformat == matrixformat::sparse ){
        
        int num_entries = -1;
        input >> rows >> columns >> num_entries;
        assert( 0 < rows && 0 < columns );
        entries.reserve( num_entries );
        
        for( int e = 0; e < num_entries; e++ ){
            
            MatrixMarketEntry mme;
            
            input >> mme.row >> mme.column >> mme.value;
            
            mme.row--; mme.column--;
            
            assert( 0 <= mme.row && 0 <= mme.column );
            assert( mme.row < rows && mme.column < columns );
            
            entries.push_back( mme );
            
        }
        
        assert( entries.size() == num_entries );
        
        if( my_matrixfeature == matrixfeature::symmetric || my_matrixfeature == matrixfeature::skew )
        {
            entries.reserve( entries.size() * 2 );
            
            double sign = my_matrixfeature == matrixfeature::symmetric ? 1. : -1;
            
            for( int e = 0; e < num_entries; e++ ){
                
                MatrixMarketEntry mme;
                mme.row = entries[e].column; 
                mme.column = entries[e].row; 
                mme.value = sign * entries[e].value;
                
                if( mme.row != mme.column ) 
                    entries.push_back( mme );
                
            }
            
        }
        
        entries.shrink_to_fit();
        
    }
    
    
    if( my_matrixformat == matrixformat::dense ){
        
        input >> rows >> columns;
        int num_entries = rows * columns;
        assert( 0 < rows && 0 < columns );
        entries.reserve( num_entries );
        
        for( int c = 0; c < columns; c++ )
        for( int r = 0; r < rows;    r++ )
        {
            MatrixMarketEntry mme;
            
            mme.row = r;
            mme.column = c;
            input >> mme.value;
            
            entries.push_back( mme );
        }
        
        assert( entries.size() == num_entries );
        
    }
        
        
    
    
    
}


void MatrixMarket::WriteSparse( std::ostream& output, const int& rows, const int& columns, const std::vector<MatrixMarket::MatrixMarketEntry>& entries, std::string comment )
{
    
    assert( rows    > 0 );
    assert( columns > 0 );
    
    output << Consts::str_magicstring << space
           << Consts::str_matrix      << space 
           << Consts::str_sparse      << space
           << Consts::str_real        << space
           << Consts::str_general     << std::endl;
           
    std::cout << '%' << comment << std::endl;
    
    output << rows << space
            << columns << space 
            << entries.size() << std::endl;
            
    for( const MatrixMarketEntry& mme : entries ){
        
        output << mme.row-1 << space
               << mme.column-1 << space
               << mme.value << std::endl;
        
    }
    
}


