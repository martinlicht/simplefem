#ifndef INCLUDEGUARD_INDEXRANGE
#define INCLUDEGUARD_INDEXRANGE


#include <cassert>

#include <vector>
#include <limits>
#include <iostream>

#include "../basic.hpp"



class IndexRange
{
    
    public:
        
        IndexRange( int, int );
        
        void check() const;
        void print( std::ostream& ) const;
        
        int min() const;
        int max() const;
        int getlength() const;
        
        bool isempty() const;
        bool contains( int ) const;
        bool contains( const IndexRange& subir ) const;
        bool operator== ( const IndexRange& ) const;
        int element2position( int ) const;
        int position2element( int ) const;
        
    private:

        int minimum;
        int maximum;
        
};

inline std::ostream& operator<<( std::ostream& os, const IndexRange& ir )
{
    ir.print( os );
    return os;
}

// static const IndexRange  NonNegativeIntegers = IndexRange( 0, std::numeric_limits<int>::max() );
// static const IndexRange  PositiveIntegers = IndexRange( 1, std::numeric_limits<int>::max() );


inline IndexRange operator|( const IndexRange& left, const IndexRange& right )
{
    return IndexRange( 
        std::min( left.min(), right.min() ),
        std::max( left.max(), right.max() )
        );
}

inline IndexRange operator&( const IndexRange& left, const IndexRange& right )
{
    return IndexRange( 
        std::max( left.min(), right.min() ), 
        std::min( left.max(), right.max() ) 
        );
}






#endif