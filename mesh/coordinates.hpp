#ifndef INCLUDEGUARD_MESH_COORDINATES_HPP
#define INCLUDEGUARD_MESH_COORDINATES_HPP


#include <array>
// #include <ostream>
#include <vector>

#include "../base/include.hpp"
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




class Coordinates final
{

    public:

        Coordinates( const Coordinates& ) = default;
        Coordinates& operator=( const Coordinates& ) = default;
        Coordinates( Coordinates&& ) noexcept = default;
        Coordinates& operator=( Coordinates&& ) noexcept = default;

        
        
        Coordinates( int dimension, int number );
        Coordinates( int dimension, int number, const std::vector<Float>& );
        ~Coordinates() noexcept = default;
        
        void check() const;
        // void print( std::ostream& ) const;
        std::string text() const;
        
        // // void lg() const { LOG << *this << nl; };
        

        // void read( std::istream& ) ;

        int getdimension() const;
        int getnumber() const;
        IndexRange getIndexRange() const;
        
        /* get/set coordinates as per point  */
        
        Float getdata( int n, int d ) const;
        void setdata( int n, int d, Float v );
        
        /* get range of coordinates */
        
        Float getmin( int d) const;
        Float getmax( int d) const;
        
        /* get/set points as vectors  */
        
        FloatVector getdata_by_vertex( int n ) const;
        FloatVector getdata_by_vertex( int n, Float s ) const;
        void setdata_by_vertex( int n, const FloatVector& input );
        void setdata_by_vertex( int n, const FloatVector& input, Float s );
        
        /* get/set coordinates as vectors  */
        
        FloatVector getdata_by_dimension( int d, Float s = 1.0 ) const;
        void setdata_by_dimension( int d, const FloatVector& value, Float s = 1.0 );
        
        /* transform all coordinates  */
        
        void scale( Float alpha );
        void scale( FloatVector alphas );
        void shift( const FloatVector& add );
        void shake_random( Float epsilon = 0.001 );
        void lineartransform( const LinearOperator& op );
        
        /* Add additional coordiantes */
        
        void append( const Coordinates& co );
        void append( const FloatVector& v );
        // void append( const std::vector<FloatVector>& );
        
        void addcapacity( int additional_capacity );
        void addcoordinates( int add_number );
        
        /* Obtain information about reference transformation of simplex */
        
        DenseMatrix getLinearPart( const IndexMap& im ) const;
        FloatVector getShiftPart( const IndexMap& im ) const;
        
        FloatVector getCenter() const;


        /* other */ 

        std::vector<Float>& raw();
        const std::vector<Float>& raw() const;
        
        std::size_t memorysize() const;
        
    private:
            
        int dimension;
        int number;
        std::vector<Float> data;

    public:

        static bool is_equal_to( const Coordinates& coords_left, const Coordinates& coords_right );

        friend inline bool operator==( const Coordinates& coords_left, const Coordinates& coords_right )
        {
            return is_equal_to( coords_left, coords_right );
        }

        friend inline bool operator!=( const Coordinates& coords_left, const Coordinates& coords_right )
        {
            return ! ( coords_left == coords_right );
        }

        template<typename Stream>
        friend inline decltype(auto) operator<<( Stream&& os, const Coordinates& co )
        {
            os << co.text();
            return std::forward<Stream>(os);
        }

// template<typename Stream>
//         friend inline std::ostream& operator<<( std::ostream& os, const Coordinates& co )
//         {
//             os << co.text(); // co.print( os );
//             return std::forward<Stream>(os);
//         }


};







#endif
