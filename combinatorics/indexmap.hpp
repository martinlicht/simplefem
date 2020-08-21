#ifndef INCLUDEGUARD_COMBINATORICS_INDEXMAP
#define INCLUDEGUARD_COMBINATORICS_INDEXMAP


#include <cassert>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <vector>

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
    
    public:
        
//         explicit IndexMap( const IndexRange& );
        IndexMap( const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const std::initializer_list<int>& );
        
//         IndexMap( const IndexRange&, const IndexRange& );
        IndexMap( const IndexRange&, const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::initializer_list<int>& );
        
        IndexMap()                              = delete;
        IndexMap( const IndexMap& )             = default;
        IndexMap& operator =( const IndexMap& ) = default;
        IndexMap( IndexMap&& )                  = default;
        IndexMap& operator =( IndexMap&& )      = default;
        virtual ~IndexMap()                     = default; // dtor virtualization is a bit shady here
        
        void check() const;
        
        void print( std::ostream&, bool embellish = true ) const;
        
        const IndexRange& getSourceRange() const;
        
        const IndexRange& getDestRange() const;
        
        
        bool isempty() const;
        
        bool isinjective() const;
        
        bool issurjective() const;
        
        bool isbijective() const;
        
        bool isstrictlyascending() const;
        
        
        // TODO
        // This interface looks like a std::vector 
        // but it should really be a mapping.
        // What's more, the return of references 
        // breaks the binding.
        
        // As a solution, the element access should only be const 
        // so it does not break the encapsulation of the class
        // and there should be explicit getter/setter methods
        
        int& at( int i ) &;
        
        const int& at( int i ) const &;
        
        int& operator[]( int i ) &;
        
        const int& operator[]( int i ) const &;
        
        const std::vector<int>& getvalues() const &;
                
        bool rangecontains( int value ) const;
        
        int preimageof( int value ) const;
        
        
        
        bool comparablewith( const IndexMap& ) const;
        
        bool equals( const IndexMap& ) const;
        
        bool less( const IndexMap& ) const;
    
    private:
        
        IndexRange src;
        IndexRange dest;
        
        std::vector<int> values;
        
};



inline IndexMap operator*( const IndexMap& leave, const IndexMap& enter )
{
    leave.check();
    enter.check();
    IndexRange src  = enter.getSourceRange();
    IndexRange dest = leave.getDestRange();

    assert( enter.getDestRange() == leave.getSourceRange() );

    IndexMap ret( src, dest, [ &leave, &enter ]( int i ) -> int { return leave[ enter[i] ]; } );

//     for( int i = src.min(); i <= src.max(); i++ )
//         ret[i] = leave[ enter[i] ];

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


inline std::ostream& operator<<( std::ostream& os, const IndexMap& im )
{
    im.check();

    im.print( os );

    return os;
}



inline IndexMap identityIndexMap( const IndexRange& ir )
{
    ir.check();

    IndexMap im( ir, ir, []( int i ) -> int { return i; } );
    im.check();

    return im;
}

inline IndexMap identityIndexMap( int low, int high )
{
    return identityIndexMap( IndexRange( low, high ) );
}



inline IndexMap expand_zero( const IndexMap& im, int p )
{
    const auto& src_range = im.getSourceRange();
    const auto& dst_range = im.getDestRange();
    
    assert( not dst_range.isempty() );
    assert( dst_range.min() == 0    );
    
    if( src_range.isempty() ) {
        
        return IndexMap( IndexRange(0,0), dst_range, {p} );
        
    } else {
        
        assert( src_range.min() == 0 );
        auto new_values = im.getvalues();
        new_values.push_back( p );
        std::sort( new_values.begin(), new_values.end() );
        const auto ret = IndexMap( IndexRange(0,src_range.max()+1), dst_range, new_values );
        assert( ret.isstrictlyascending() );
        return ret;
        
    }
}


inline IndexMap expand_one( const IndexMap& im, int p )
{
    const auto& src_range = im.getSourceRange();
    const auto& dst_range = im.getDestRange();
    
    assert( not dst_range.isempty() );
    assert( dst_range.min() == 0    );
    
    if( src_range.isempty() ) {
        
        return IndexMap( IndexRange(1,1), dst_range, {p} );
        
    } else {
        
        assert( src_range.min() == 1 );
        auto new_values = im.getvalues();
        new_values.push_back( p );
        std::sort( new_values.begin(), new_values.end() );
        const auto ret = IndexMap( IndexRange(1,src_range.max()+1), dst_range, new_values );
        assert( ret.isstrictlyascending() );
        return ret;
        
    }
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
//             unreachable();
// }



#endif
