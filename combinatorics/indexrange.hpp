#ifndef INCLUDEGUARD_COMBINATORICS_INDEXRANGE
#define INCLUDEGUARD_COMBINATORICS_INDEXRANGE



#include <vector>
#include <limits>
#include <iostream>

#include "../basic.hpp"

/*****
**
**  This class models a range of indices, i.e. integers,
**  of the form \{ min, min+1, ..., max \}.
**  The index range may be empty.
**
******/


class IndexRange
{
    
    public:
        
        IndexRange( int from, int to );
        
        void check() const;
        void print( std::ostream& ) const;
        
        int min() const;
        int max() const;
        int cardinality() const;
        int getlength() const; /* synonym to cardinality */
        
        bool isempty() const;
        bool contains( int element ) const;
        bool contains( const IndexRange& subir ) const;
        bool compare( const IndexRange&  ) const;
        int element2position( int element ) const;
        int position2element( int position ) const;
        
        /* For each semantics */ 
        
        class IndexRangeIterator {
            
            private: 
            
                int value;
                
            public:
                
                IndexRangeIterator(int value) : value(value) {};
                inline int operator*() const { return value; };
                inline IndexRangeIterator& operator++() { ++value; return *this; };
                inline IndexRangeIterator operator++( int ) { return (++value)-1; };
                inline bool operator==( const IndexRangeIterator& irit ) const 
                    { return value == irit.value; };
                inline bool operator!=( const IndexRangeIterator& irit ) const 
                    { return value != irit.value; };
                    
        };
        
        inline const IndexRangeIterator begin() const { return IndexRangeIterator(minimum); };
        inline const IndexRangeIterator end() const { return IndexRangeIterator(maximum+1); };
        
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


// inline IndexRange operator|( const IndexRange& left, const IndexRange& right )
// {
//     return IndexRange( 
//         std::min( left.min(), right.min() ),
//         std::max( left.max(), right.max() )
//         );
// }
// 
// inline IndexRange operator&( const IndexRange& left, const IndexRange& right )
// {
//     return IndexRange( 
//         std::max( left.min(), right.min() ), 
//         std::min( left.max(), right.max() ) 
//         );
// }

inline bool operator== ( const IndexRange& ir1, const IndexRange& ir2 )
{
    return ir1.compare( ir2 );
}

inline bool operator!= ( const IndexRange& ir1, const IndexRange& ir2 )
{
    return !( ir1 == ir2 );
}




#endif