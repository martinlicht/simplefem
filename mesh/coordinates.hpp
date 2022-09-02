#ifndef INCLUDEGUARD_MESH_COORDINATES_HPP
#define INCLUDEGUARD_MESH_COORDINATES_HPP


// #include <cassert>
#include <ostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../dense/densematrix.hpp"



/*******************
****  
****  Class for Coordiante Collections 
****  
****  - Most basic functionality
****  
****  
****  
****  
*******************/




class Coordinates
{

    public:

        Coordinates( const Coordinates& ) = default;
        Coordinates& operator=( const Coordinates& ) = default;
        Coordinates( Coordinates&& ) = default;
        Coordinates& operator=( Coordinates&& ) = default;

        
        
        Coordinates( int dimension, int number );
        Coordinates( int dimension, int number, const std::vector<Float>& );
        virtual ~Coordinates();
        
        void check() const;
        void print( std::ostream& ) const;
        std::string text() const;
        
        // // void lg() const { LOG << *this << nl; };
        

        void read( std::istream& ) ;

        int getdimension() const;
        int getnumber() const;
        IndexRange getIndexRange() const;
        
        /* get/set coordinates as per point  */
        
        Float getdata( int n, int d) const;
        void setdata( int n, int d, Float v );
        
        /* get/set points as vectors  */
        
        FloatVector getvectorclone( int n ) const;
        FloatVector getvectorclone( int n, Float s ) const;
        void loadvector( int n, const FloatVector& value );
        void loadvector( int n, const FloatVector& value, Float s );
        
        /* get/set coordinates as vectors  */
        
        FloatVector getdimensionclone( int d, Float s = 1.0 ) const;
        void loaddimension( int d, const FloatVector& value, Float s = 1.0 );
        
        /* transform all coordinates  */
        
        void scale( Float );
        void shift( const FloatVector& );
        void lineartransform( const LinearOperator& );
        
        /* Add additional coordiantes */
        
        void append( const Coordinates& );
        void append( const FloatVector& );
        void append( const std::vector<FloatVector>& );
        
        void addcapacity( int );
        void addcoordinates( int );
        
        /* Obtain information about reference transformation of simplex */
        
        DenseMatrix getLinearPart( const IndexMap& ) const;
        FloatVector getShiftPart( const IndexMap& ) const;
        
        FloatVector getCenter() const;
        
    private:
            
        int dimension;
        int number;
        std::vector<Float> data;
        
};

bool compare( const Coordinates& coords_left, const Coordinates& coords_right, Float tolerance = 0.00001 );

inline std::ostream& operator<<( std::ostream& os, const Coordinates& co )
{
    co.print( os );
    return os;
}

inline bool operator==( const Coordinates& coords_left, const Coordinates& coords_right )
{
    return compare( coords_left, coords_right, 0.00001 );
}

inline bool operator!=( const Coordinates& coords_left, const Coordinates& coords_right )
{
    return ! ( coords_left == coords_right );
}


#endif
