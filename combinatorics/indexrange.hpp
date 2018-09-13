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
        void print( std::ostream&, bool embellish = false ) const;
        
        int min() const;
        int max() const;
        int cardinality() const;
        
        bool isempty() const;
        bool contains( int element ) const;
        bool contains( const IndexRange& subir ) const;
        bool operator== ( const IndexRange&  ) const;
        int element2position( int element ) const;
        int position2element( int position ) const;
        
        /* For each semantics */ 
        
        class IndexRangeIterator {
            
            private: 
            
                int value;
                
            public:
                
                explicit IndexRangeIterator(int value) : value(value) 
                { }
                
                inline int operator*() const
                {
                    return value;                
                }
                
                inline IndexRangeIterator& operator++()
                {
                    ++value; 
                    return *this;
                } // pre-increment
                
                inline IndexRangeIterator operator++( int )
                {
                    return IndexRangeIterator( value++ );   
                } // post-increment
                
                inline bool operator==( const IndexRangeIterator& irit ) const 
                { 
                    return value == irit.value;
                }
                
                inline bool operator!=( const IndexRangeIterator& irit ) const 
                { 
                    return value != irit.value;
                }
                    
        };
        
        inline const IndexRangeIterator begin() const
        {
            return IndexRangeIterator(minimum);
        }
        
        inline const IndexRangeIterator end() const
        {
            return IndexRangeIterator(maximum+1);
        }
        
        /* enum class */
        
        // TODO: make general conceptions about how to handle the output of objects 
        // in these modules. Generally, we would like to control the output format by some parameter 
        // given to each print function. We can assume that the parameters belong to some enum class 
        // defined within a class declaration and are specifcally tailored to each class. 
        // They are a purely optional argument for the print method and may be skipped at convenience.
        
    private:

        int minimum;
        int maximum;
        
};

inline std::ostream& operator<<( std::ostream& os, const IndexRange& ir )
{
    ir.print( os );
    return os;
}

// static const IndexRange  NonNegativeIntegers = IndexRange( 0, std::numeric_limits<int>::max()-1 );
// 
// static const IndexRange  PositiveIntegers = IndexRange( 1, std::numeric_limits<int>::max()-1 );
// 
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


// inline bool operator==( const IndexRange& left, const IndexRange& right )
// {
//     if( left.isempty() )
//         return left.isempty() == right.isempty();
//     else 
//         return left.min() == right.min() && left.max() == right.max();
// }



#endif
