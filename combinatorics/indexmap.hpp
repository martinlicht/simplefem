#ifndef INCLUDEGUARD_COMBINATORICS_INDEXMAP_HPP
#define INCLUDEGUARD_COMBINATORICS_INDEXMAP_HPP


#include <cassert>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <string>
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
        
        /* Constructors */
        
        IndexMap( const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const std::initializer_list<int>& );
        
        IndexMap( const IndexRange&, const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::initializer_list<int>& );
        
        /* standard interface */ 
        
        IndexMap()                              = delete;
        IndexMap( const IndexMap& )             = default;
        IndexMap& operator =( const IndexMap& ) = default;
        IndexMap( IndexMap&& )                  = default;
        IndexMap& operator =( IndexMap&& )      = default;
        virtual ~IndexMap()                     = default; // dtor virtualization is a bit shady here
        
        /* standard methods */
        
        void check() const;
        
        std::string text( bool embellish = true ) const;
        
        void print( std::ostream&, bool embellish = true ) const;

        void lg() { LOG << *this << nl; };

        
        /* OTHER METHODS */
        
        const IndexRange& getSourceRange() const;
        
        const IndexRange& getTargetRange() const;
        
        bool isempty() const;
        
        bool isinjective() const;
        
        bool issurjective() const;
        
        bool isbijective() const;
        
        bool isstrictlyascending() const;
        
        
        // TODO
        // This interface looks like a std::vector 
        // but it should really be a mapping.
        // What's more, the return of references 
        // may break the codomain, so to speak.
        
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
        
        
        
        // IndexMap skip( int i ) const;
        // IndexMap shiftup() const;



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
    IndexRange dest = leave.getTargetRange();

    assert( enter.getTargetRange() == leave.getSourceRange() );

    IndexMap ret( src, dest, [ &leave, &enter ]( int i ) -> int { return leave[ enter[i] ]; } );

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



// inline IndexMap toss_entry( const IndexMap& original, int p )
// {
    
//     assert( not original.getSourceRange().isempty() );
//     assert( original.getSourceRange().contains(p) );
    
//     IndexRange new_source_range = IndexRange( original.getSourceRange().min()+1, original.getSourceRange().max() );

//     IndexMap im( new_source_range, original.getTargetRange(), [p,original]( int i ) -> int { 
//         if( i <= p ) 
//             return original.at(p-1);
//         else 
//             return original.at(p);
//     } );

//     // IndexMap im( new_source_range, original.getTargetRange(), 0 );
//     // for( int j = 1; j <= p; j++ ) im[j] = original[j-1];
//     // for( int j = p+1; j < original.getSourceRange().max(); j++ ) im[j] = original[j];
    
//     return im;
// }





IndexMap expand_zero( const IndexMap& im, int p );

IndexMap expand_one( const IndexMap& im, int p );



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
