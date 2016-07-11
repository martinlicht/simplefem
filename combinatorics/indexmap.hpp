#ifndef INCLUDEGUARD_INDEXMAP
#define INCLUDEGUARD_INDEXMAP


#include <cassert>

#include <vector>
#include <limits>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"






class IndexMap
{
    
    private:
        
        IndexRange src;
        IndexRange dest;
        
        std::vector<int> values;
        
    public:
        
        IndexMap( const IndexRange& );
        IndexMap( const IndexRange&, const IndexRange& );
        IndexMap( const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::function<int(int)>& );
        
        void check() const;
        void print( std::ostream& ) const;
        
        const IndexRange& getSourceRange() const;
        const IndexRange& getDestRange() const;
        
        bool isinjective() const;
        bool issurjective() const;
        bool isbijective() const;
        bool isstrictlyascending() const;
        
        int& operator[]( int i );
        const int& operator[]( int i ) const;
        const std::vector<int>& getvalues() const;
        
        IndexMap inverse() const;
        
        IndexMap skip( int ) const;
        IndexMap attachbefore( int ) const;
        
        bool comparablewith( const IndexMap& ) const;
        bool equals( const IndexMap& ) const;
        bool less( const IndexMap& ) const;
        
};



inline IndexMap operator*( const IndexMap& leave, const IndexMap& enter )
{
    IndexRange src  = enter.getSourceRange();
    IndexRange dest = leave.getDestRange();

    assert( enter.getDestRange() == leave.getSourceRange() );

    IndexMap ret( src, dest );

    for( int i = src.min(); i <= src.max(); i++ )
        ret[i] = leave[ enter[i] ];

    return ret;
}

inline bool operator==( const IndexMap& left, const IndexMap& right )
{
    assert( left.comparablewith( right ) );
    return left.equals( right );
}

inline bool operator!=( const IndexMap& left, const IndexMap& right )
{
    assert( left.comparablewith( right ) );
    return !( left.equals( right ) );
}

inline bool operator<( const IndexMap& left, const IndexMap& right )
{
    assert( left.comparablewith( right ) );
    return left.less( right );
}


inline std::ostream& operator<<( std::ostream& os, IndexMap im )
{
    im.print( os );
    return os;
}



inline IndexMap identityIndexMap( int low, int high )
{
    IndexRange ir( low, high );
    IndexMap im( ir, ir );
    for( int i = low; i <= high; i++ )
            im[i] = i;
    return im;
}






#endif