#ifndef INCLUDEGUARD_COMBINATORICS_PERMUTATION
#define INCLUDEGUARD_COMBINATORICS_PERMUTATION


#include <vector>
#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


/********************
**** 
****  Class that describes permutations over an IndexRange 
****  defines the sign of a permutation
**** 
********************/

class Permutation
: public IndexMap
{

        public:
        
                explicit Permutation( const IndexRange& ir );
                Permutation( const IndexRange& ir, const std::vector<int>& );
                
                void check() const;

                void print( std::ostream&, bool embellish = false ) const;
                
                IndexRange getIndexRange() const;

                const std::vector<int>& getvalues() const;
                
                bool getsign();
        
        private:
    
                int sign;
};




inline std::ostream& operator<<( std::ostream& os, const Permutation& mi )
{
        mi.check();
        mi.print( os );
        return os;
}

inline Permutation& operator+=( Permutation& left, int right )
{
        left.check();
        left.add( right );
        return left;
}

inline Permutation& operator-=( Permutation& left, int right )
{
        left.check();
        left.sub( right );
        return left;
}

inline Permutation& operator+=( Permutation& left, const Permutation& right )
{
        left.check();
        right.check();
        left.add( right );
        return left;
}

inline Permutation& operator-=( Permutation& left, const Permutation& right )
{
        left.check();
        right.check();
        left.sub( right );
        return left;
}


inline Permutation operator+( const Permutation& left, int right )
{
        left.check();
        Permutation ret = left;
        ret += right;
        ret.check();
        return ret;
}

inline Permutation operator-( const Permutation& left, int right )
{
        left.check();
        Permutation ret = left;
        ret -= right;
        ret.check();
        return ret;
}

inline Permutation operator+( const Permutation& left, const Permutation& right )
{
        left.check();
        right.check();
        Permutation ret = left;
        ret += right;
        ret.check();
        return ret;
}

inline Permutation operator-( const Permutation& left, const Permutation& right )
{
        left.check();
        right.check();
        Permutation ret = left;
        ret -= right;
        ret.check();
        return ret;
}


inline bool operator==( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return it.equals( mi );
}
                
inline bool operator!=( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return ! ( it == mi );
}

inline bool operator<( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return it.less( mi );
}
                
inline bool operator>( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return mi < it;
}
                
inline bool operator<=( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return it < mi || it == mi;
}
                
inline bool operator>=( const Permutation& it, const Permutation& mi)
{
        it.check();
        mi.check();
        return it > mi || it == mi;
}
        
        
inline int absolute( const Permutation& it )
{
        return it.absolute();
}


inline int factorial( const Permutation& it )
{
        return it.factorial();
}




                
#endif
