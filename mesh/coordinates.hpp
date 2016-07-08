#ifndef INCLUDEGUARD_COORDINATES
#define INCLUDEGUARD_COORDINATES


#include <vector>
#include <iostream>
#include <cassert>
#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
// #include "../.hpp"
// #include "../.hpp"



class Coordinates
{

    public:

        Coordinates( int, int );

        void check() const;
        void print( std::ostream& ) const;

        void read( std::istream& ) ;

        int getdimension() const;
        int getnumber() const;
        
        Float getdata( int, int ) const;
        void setdata( int, int, Float );
        
        FloatVector getvector( int ) const;
        void setvector( int, const FloatVector& );
        
        void scale( Float );
        void shift( const FloatVector& );
        void lineartransform( const LinearOperator& );
        
        void append( const Coordinates& );
        void append( const FloatVector& );
            
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