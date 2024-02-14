#ifndef INCLUDEGUARD_OPERATOR_FLOATVECTOR_HPP
#define INCLUDEGUARD_OPERATOR_FLOATVECTOR_HPP

#include <functional>
#include <initializer_list>
// #include <ostream>
#include <vector>

#include "../basic.hpp"

class LinearOperator;

/**********************
***
***  Describes a vector of Floating point numbers 
***  offers basic arithmetic operations
***
**********************/


class FloatVector
{
    
    public:
        
        /* Constructors */
        
        explicit FloatVector( int dim, Float initivalue = notanumber );
        
        explicit FloatVector( const FloatVector&, Float scaling );
        explicit FloatVector( FloatVector&&, Float scaling );
        
        explicit FloatVector( const std::vector<Float>&, Float scaling = 1. );
        explicit FloatVector( const std::vector<int>&, Float scaling = 1. );

                 FloatVector( const std::initializer_list<Float>& l );
        
        explicit FloatVector( int dimension, const std::function<Float(int)>& generator, Float scaling = 1. );


        /* standard interface */ 
        
        FloatVector() = delete;
        FloatVector( const FloatVector& );
        FloatVector( FloatVector&& );
        FloatVector& operator=( const FloatVector& vec );
        FloatVector& operator=( FloatVector&& vec );

        virtual ~FloatVector();
        
        
        
        /* standard methods */

        void check() const;
        
        std::string text() const; 

        std::string data_as_text( bool indexed = true, bool rowwise = false ) const; 

        // void print( std::ostream& ) const;
        
        // void lg() const { LOG << text() << nl; };
        
        
        /* OTHER METHODS */

        /* Cloning */

        FloatVector clone() const;


        /* information and data access */
        
        HOTCALL int getdimension() const;
        
        Float setentry( int, Float );
        
        Float getentry( int ) const;
        

        Float& at( int ) &;
        
        const Float& at( int ) const &;
        
        Float& operator[]( int ) &;
        
        const Float& operator[]( int ) const &;
        
        const std::vector<Float> getdata() const;
        
        
        /* load values */
        
        void setentries( Float );
        
        void random();
        
        void random_within_range( Float min, Float max );
        
        void to_absolute();
        
        void zero();
        
        void clear();
        
        void clear_if( const std::vector<bool>& mask );
        
        void clear_unless( const std::vector<bool>& mask );
        
        
        /* scale and shift */
        
        FloatVector& normalize();
        
        FloatVector& normalize( const LinearOperator& );
        
        FloatVector& scale( Float );
        
        FloatVector& scaleinverse( Float );
        
        FloatVector& shift( Float );
        
        FloatVector& shiftnegative( Float );
        
        
        /* slices */
        
        FloatVector getslice( int, int ) const;
        
        void setslice( int, int, Float );
        
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
        
        Float scalarproductwith( const FloatVector&, const std::vector<bool>& ) const;
        
        
        /* Calculations */
        
        Float min() const;
        
        Float average() const;
        
        Float sum() const;
        
        Float norm() const;
        
        Float norm_sq() const;
        
        Float norm( const LinearOperator& ) const;
        
        Float norm_sq( const LinearOperator& ) const;
        
        Float maxnorm() const;
        
        Float sumnorm() const;
        
        Float l2norm() const;
        
        Float lpnorm( Float p, Float inner_weight = 1. ) const;
        
        
        
        
        /* Investigations */
        
        bool isfinite() const;
        
        bool iszero() const;
        
        bool ispositive() const;
        
        bool isnegative() const;
        
        bool isnonnegative() const;
        
        bool isnonpositive() const;
        
        
        bool is_numerically_small( Float threshold = desired_closeness ) const;
        


        
        
        /* Raw access */
        
        Float* raw();
        
        const Float* raw() const;
        
        
        
        
        /* Memory size */
        
        std::size_t memorysize() const;
        
        
        
        
        /* For each semantics */ 
        
        class ConstIterator {
            
            friend FloatVector;
            
            private: 
            
                const FloatVector* parentvector;
                int index;
                
                explicit ConstIterator( const FloatVector* parentvector, int index) : parentvector(parentvector), index(index) 
                { }
                
            public:
                
                inline Float operator*() const
                {
                    return parentvector->at(index);
                }
                
                inline ConstIterator operator++()
                {
                    ++index; 
                    return *this;
                } // pre-increment
                
                inline ConstIterator operator++( int )
                {
                    return ConstIterator( parentvector, index++ );   
                } // post-increment
                
                inline bool operator!=( const ConstIterator& other ) const 
                { 
                    assert( other.parentvector == parentvector );
                    return index != other.index;
                }
                    
        };
        
        class Iterator {
            
            friend FloatVector;
            
            private: 
            
                FloatVector* parentvector;
                int index;
                
                explicit Iterator( FloatVector* parentvector, int index) : parentvector(parentvector), index(index) 
                { }
                
            public:
                
                inline Float& operator*() const
                {
                    return parentvector->at(index);
                }
                
                inline Iterator operator++()
                {
                    ++index; 
                    return *this;
                } // pre-increment
                
                inline Iterator operator++( int )
                {
                    return Iterator( parentvector, index++ );   
                } // post-increment
                
                inline bool operator!=( const Iterator& other ) const 
                { 
                    assert( other.parentvector == parentvector );
                    return index != other.index;
                }
                    
        };
        
        
        
        inline ConstIterator begin() const
        {
            return ConstIterator(this,0);
        }
        
        inline ConstIterator end() const
        {
            return ConstIterator(this,dimension);
        }

        inline Iterator begin()
        {
            return Iterator(this,0);
        }
        
        inline Iterator end()
        {
            return Iterator(this,dimension);
        }

        
        
        
    private:

        int dimension;
        Float* pointer;

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
// inline std::ostream& operator<<( std::ostream& out, const FloatVector& vec )
template<typename Stream>
// inline std::ostream& operator<<( std::ostream& out, const FloatVector& vec )
inline Stream& operator<<( Stream&& out, const FloatVector& vec )
{
    out << vec.text(); // vec.print( out );
    return out;
}




inline FloatVector unitvector( int d, int i )
{
    FloatVector ret( d, 0.0 );
    ret[ i ] = 1.0;
    return ret;
}




// Base function to end recursion
inline FloatVector concatFloatVector( const FloatVector& vec ) {
    return vec;
}

// Template function to concatenate strings
template <typename... Args>
inline FloatVector concatFloatVector( const FloatVector& first, Args... args) {
    auto others = concatFloatVector(args...);
    FloatVector ret( first.getdimension() + others.getdimension() );
    ret.setslice( 0, first );
    ret.setslice( first.getdimension(), others );
    return ret;
}

//
FloatVector interlace( const FloatVector& first, const FloatVector& second );



#endif
