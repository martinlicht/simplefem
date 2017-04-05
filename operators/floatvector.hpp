#ifndef INCLUDEGUARD_FLOATVECTOR
#define INCLUDEGUARD_FLOATVECTOR



// #include <cassert>

#include <ostream>
#include <vector>
#include <functional>

#include "../basic.hpp"


/**********************
***
***  Describes a vector of Floating point numbers 
***  offers basic arithmetic operations
***
**********************/


class FloatVector
{

    public:

        explicit FloatVector( int dim );
        FloatVector( int dim, Float initivalue );
        FloatVector( const FloatVector& );
        explicit FloatVector( const std::vector<Float>& );
        FloatVector( const std::vector<Float>&, Float );
        FloatVector( const FloatVector&, Float scaling );
        FloatVector( int dimension, const std::function<Float(int)>& generator );
        FloatVector( int dimension, const std::function<Float(int)>& generator, Float scaling );
        // virtual ~FloatVector();

        void check() const;

        void print( std::ostream& ) const;
        void printplain( std::ostream& ) const;
        
        /* information and data access */
        
        int getdimension() const;
        
        Float setentry( int, Float );
        Float getentry( int ) const;

        Float& at( int );
        const Float& at( int ) const;
        Float& operator[]( int );
        const Float& operator[]( int ) const;

        const std::vector<Float>& getdata() const;
        
        /* basic manipulation */
        
        void zero();
        void random();
        void scale( Float );

        /* slices */
        
        FloatVector getslice( int, int ) const;
        void setslice( int, const FloatVector& );
        void addslice( int, const FloatVector&, Float );
        
        /* arithmetics and assignments */
        
        void copydatafrom( const FloatVector& );
        void copydatafrom( Float, const FloatVector& );

        void generatedatafrom( const std::function<Float(int)>& );
        void generatedatafrom( Float, const std::function<Float(int)>& );

        void adddatafrom( const FloatVector& );
        void adddatafrom( Float, const FloatVector& );
        void adddatafrom( Float, Float, const FloatVector& );

        Float scalarproductwith( const FloatVector& ) const;
        
        /* Calculations */
        
        Float norm() const;
        Float maxnorm() const;
        Float lpnorm( Float ) const;
        
        
    private:

        // int dimension;
        std::vector<Float> data;

};



/* 
 *
 * Overload Operators  
 *
 */

inline FloatVector operator+( const FloatVector& vec )
{
    return vec;
}

inline FloatVector operator-( const FloatVector& vec )
{
    FloatVector ret( vec, -1. ); 
    return ret;
}

/* Scaling operations */ 

inline FloatVector operator*=( FloatVector& vec, Float scaling )
{
    vec.scale( scaling );
    return vec;
}

inline FloatVector operator/=( FloatVector& vec, Float scaling )
{
    vec.scale( 1. / scaling );
    return vec;
}

/* Adding and Subtracting operations */ 

inline FloatVector operator+=( FloatVector& vec, const FloatVector& src )
{
    vec.adddatafrom( src );
    return vec;
}

inline FloatVector operator-=( FloatVector& vec, const FloatVector& src )
{
    vec.adddatafrom( -1., src );
    return vec;
}

inline FloatVector operator+( const FloatVector& left, const FloatVector& right )
{
    FloatVector vec( left );
    vec += right;
    return vec;
}

inline FloatVector operator-( const FloatVector& left, const FloatVector& right )
{
    FloatVector vec( left );
    vec -= right; 
    return vec;
}

inline FloatVector operator*( Float s, const FloatVector& vec )
{
    FloatVector ret(vec, s );
    return ret;
}

inline FloatVector operator/( const FloatVector& vec, Float s )
{
    FloatVector ret(vec, 1. / s );
    return ret;
}

inline Float operator*( const FloatVector& left, const FloatVector& right )
{
    // assert( left.getdimension() == right.getdimension() );
    // Float ret = 0.;
    // for( int p = 0; p < left.getdimension(); p++ )
        // ret += left.getentry(p) * right.getentry(p);
    // return ret;
    return left.scalarproductwith( right );
}


/* Output stream notation */
inline std::ostream& operator<<( std::ostream& out, const FloatVector& vec )
{
    vec.print( out );
    return out;
}
 


#endif