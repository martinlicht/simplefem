#ifndef INCLUDEGUARD_COMBINATORICS_MULTIINDEX_HPP
#define INCLUDEGUARD_COMBINATORICS_MULTIINDEX_HPP


#include <functional>
#include <initializer_list>
#include <iostream>
#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "indexmap.hpp"


/********************
**** 
****  Class that describes multiindices over an IndexRange 
****  basic arithmetic functionality
**** 
********************/

class MultiIndex final
: public IndexMap
{

        public:
        
                explicit MultiIndex( const IndexRange& ir );
                MultiIndex( const IndexRange& ir, const std::vector<int>& );
                MultiIndex( const IndexRange&, const std::function<int(int)>& );
                MultiIndex( const IndexRange&, const std::initializer_list<int>& );
                MultiIndex( const MultiIndex& )             = default;
                MultiIndex& operator =( const MultiIndex& ) = default;
                MultiIndex( MultiIndex&& )                  = default;
                MultiIndex& operator =( MultiIndex&& )      = default;
                virtual ~MultiIndex() = default;
                
                void check() const;

                void print( std::ostream&, bool embellish = false ) const;
                
                IndexRange getIndexRange() const;

                const std::vector<int>& getvalues() const;
        
                
//                 const int& at(int) const;
// 
//                 int& at(int);
// 
//                 const int& operator[](int) const;
// 
//                 int& operator[](int);
                

                int absolute() const;

                int factorial() const;
                
                int min() const;

                int max() const;
                

                void add( int );

                void sub( int );
                
                void add( int, int );

                void sub( int, int );
                
                void add( const MultiIndex& );

                void sub( const MultiIndex& );
                

                bool comparablewith( const MultiIndex& ) const;

                bool less( const MultiIndex& ) const;

                bool equals( const MultiIndex& ) const;
                
//         private:
// 
//                 IndexRange range;
//                 
//                 std::vector<int> values;
                
};




inline std::ostream& operator<<( std::ostream& os, const MultiIndex& mi )
{
        mi.check();
        mi.print( os );
        return os;
}

inline MultiIndex& operator+=( MultiIndex& left, int right )
{
        left.check();
        left.add( right );
        return left;
}

inline MultiIndex& operator-=( MultiIndex& left, int right )
{
        left.check();
        left.sub( right );
        return left;
}

inline MultiIndex& operator+=( MultiIndex& left, const MultiIndex& right )
{
        left.check();
        right.check();
        left.add( right );
        return left;
}

inline MultiIndex& operator-=( MultiIndex& left, const MultiIndex& right )
{
        left.check();
        right.check();
        left.sub( right );
        return left;
}


inline MultiIndex operator+( const MultiIndex& left, int right )
{
        left.check();
        MultiIndex ret = left;
        ret += right;
        ret.check();
        return ret;
}

inline MultiIndex operator-( const MultiIndex& left, int right )
{
        left.check();
        MultiIndex ret = left;
        ret -= right;
        ret.check();
        return ret;
}

inline MultiIndex operator+( const MultiIndex& left, const MultiIndex& right )
{
        left.check();
        right.check();
        MultiIndex ret = left;
        ret += right;
        ret.check();
        return ret;
}

inline MultiIndex operator-( const MultiIndex& left, const MultiIndex& right )
{
        left.check();
        right.check();
        MultiIndex ret = left;
        ret -= right;
        ret.check();
        return ret;
}


inline bool operator==( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return it.equals( mi );
}
                
inline bool operator!=( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return ! ( it == mi );
}

inline bool operator<( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return it.less( mi );
}
                
inline bool operator>( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return mi < it;
}
                
inline bool operator<=( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return it < mi || it == mi;
}
                
inline bool operator>=( const MultiIndex& it, const MultiIndex& mi)
{
        it.check();
        mi.check();
        return it > mi || it == mi;
}
        
        
inline int absolute( const MultiIndex& it )
{
        return it.absolute();
}


inline int factorial( const MultiIndex& it )
{
        return it.factorial();
}




                
#endif
