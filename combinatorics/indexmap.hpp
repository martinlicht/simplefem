#ifndef INCLUDEGUARD_INDEXMAP
#define INCLUDEGUARD_INDEXMAP

#include <vector>
#include <limits>
#include <iostream>

#include "../basic.hpp"
#include "indexrange.hpp"


/*****
**
**  This class models a mapping of indices.
**  It has a 'source' and 'target' range.
**  The mapping may be empty.
**
******/



class IndexMap
{
    
    private:
        
        IndexRange src;
        IndexRange dest;
        
        std::vector<int> values;
        
    public:
        
        explicit IndexMap( const IndexRange& );
        IndexMap( const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const std::initializer_list<int>& );
        
        IndexMap( const IndexRange&, const IndexRange& );
        IndexMap( const IndexRange&, const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::initializer_list<int>& );
        
        void check() const;
        void print( std::ostream& ) const;
        
        const IndexRange& getSourceRange() const;
        const IndexRange& getDestRange() const;
        
        bool isempty() const;
        bool isinjective() const;
        bool issurjective() const;
        bool isbijective() const;
        bool isstrictlyascending() const;
        
        int& at( int i );
        const int& at( int i ) const;
        int& operator[]( int i );
        const int& operator[]( int i ) const;
        const std::vector<int>& getvalues() const;
        
//         IndexMap inverse() const;
        
//         IndexMap skip( int ) const;
//         IndexMap attachbefore( int ) const;
        
        bool rangecontains( int ) const;
        int rangeposition( int ) const;
        
        bool comparablewith( const IndexMap& ) const;
        bool equals( const IndexMap& ) const;
        bool less( const IndexMap& ) const;
        
};



inline IndexMap operator*( const IndexMap& leave, const IndexMap& enter )
{
    leave.check();
    enter.check();
    IndexRange src  = enter.getSourceRange();
    IndexRange dest = leave.getDestRange();

    assert( enter.getDestRange() == leave.getSourceRange() );

    IndexMap ret( src, dest );

    for( int i = src.min(); i <= src.max(); i++ )
        ret[i] = leave[ enter[i] ];

    ret.check();
    return ret;
}

inline bool operator==( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.comparablewith( right ) );
    return left.equals( right );
}

inline bool operator!=( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.comparablewith( right ) );
    return !( left.equals( right ) );
}

inline bool operator<( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.comparablewith( right ) );
    return left.less( right );
}


inline std::ostream& operator<<( std::ostream& os, IndexMap im )
{
    im.check();
    im.print( os );
    return os;
}



inline IndexMap identityIndexMap( const IndexRange& ir )
{
    ir.check();
    IndexMap im( ir, ir );
    for( int i = ir.min(); i <= ir.max(); i++ )
        im[i] = i;
    im.check();
    return im;
}

inline IndexMap identityIndexMap( int low, int high )
{
    return identityIndexMap( IndexRange( low, high ) );
}


// inline int fehlstelle( const IndexMap& sub, const IndexMap& super )
// {
//     sub.check();
//     super.check();
//     assert( sub.getSourceRange().getlength() == super.getSourceRange().getlength() + 1 );
//     for( int i : sub.getSourceRange() )
//         assert( super.rangecontains(i) );
//     for( int j : super.getSourceRange() )
//         if( ! sub.rangecontains(j) )
//             return j;
//         else 
//             assert(false);
// }



#endif