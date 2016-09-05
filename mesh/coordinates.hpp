#ifndef INCLUDEGUARD_COORDINATES
#define INCLUDEGUARD_COORDINATES


#include <vector>
#include <iostream>
// #include <cassert>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/densematrix.hpp"
#include "../operators/linearoperator.hpp"
// #include "../.hpp"
// #include "../.hpp"



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

        Coordinates( int, int );

        void check() const;
        void print( std::ostream& ) const;

        void read( std::istream& ) ;

        int getdimension() const;
        int getnumber() const;
        IndexRange getIndexRange() const;
        
        /* get/set coordinates as per entry  */
        
        Float getdata( int, int ) const;
        void setdata( int, int, Float );
        
        /* get/set coordinates as vectors  */
        
        FloatVector getvectorclone( int ) const;
        FloatVector getvectorclone( int, Float ) const;
        void loadvector( int, const FloatVector& );
        void loadvector( int, const FloatVector&, Float );
        
        /* transform all coordinates  */
        
        void scale( Float );
        void shift( const FloatVector& );
        void lineartransform( const LinearOperator& );
        
        /* Add additional coordiantes */
        
        void append( const Coordinates& );
        void append( const FloatVector& );

        /* Obtain information about reference transformation of simplex */
        
        DenseMatrix getLinearPart( const IndexMap& ) const;
        FloatVector getShiftPart( const IndexMap& ) const;
        
    private:
            
        int dimension;
        int number;
        std::vector<Float> data;
        
};

inline std::ostream& operator<<( std::ostream& os, const Coordinates& ir )
{
    ir.print( os );
    return os;
}



#endif