#ifndef INCLUDEGUARD_FLOATVECTOR_HPP
#define INCLUDEGUARD_FLOATVECTOR_HPP



#include <ostream>
#include <vector>
#include <cassert>
#include "../basic.hpp"


class FloatVector
{
	
	public:

		explicit FloatVector( int dim );
		FloatVector( const FloatVector& );
		explicit FloatVector( const FloatVector&, Float scaling );
		// virtual ~FloatVector();
		
		void check() const;
		
		void print( std::ostream& ) const;
		
		void zero();
		void random();
		void scale( Float );
		void copydatafrom( const FloatVector&, Float scaling = 1. );
		void adddatafrom( const FloatVector&, Float scaling = 1. );
		
		Float setentry( int, Float );
		Float getentry( int ) const;
		
		int getdimension() const;
	
	private:
	
		int dimension;
		std::vector<Float> data;
	
};



/* 
 *
 * Overload Operators  
 *
 */
	
inline FloatVector operator+( const FloatVector& V )
{
	return V;
}

inline FloatVector operator-( const FloatVector& V )
{
	FloatVector ret( V, -1. ); 
	return ret;
}

/* Scaling operations */ 

inline FloatVector operator*=( FloatVector& V, Float scaling )
{
	V.scale( scaling );
	return V;
}
	
inline FloatVector operator/=( FloatVector& V, Float scaling )
{
	V.scale( 1. / scaling );
	return V;
}

/* Adding and Subtracting operations */ 

inline FloatVector operator+=( FloatVector& V, const FloatVector& Vsrc )
{
	V.adddatafrom( Vsrc );
	return V;
}

inline FloatVector operator-=( FloatVector& V, const FloatVector& Vsrc )
{
	V.adddatafrom( Vsrc, -1 );
	return V;
}	
	
inline FloatVector operator+( const FloatVector& V1, const FloatVector& V2 )
{
	FloatVector V( V1 );
	V += V2;
	return V;
}

inline FloatVector operator-( const FloatVector& V1, const FloatVector& V2 )
{
	FloatVector V( V1 );
	V -= V2; 
	return V;
}	

inline FloatVector operator*( Float s, const FloatVector& vec )
{
	FloatVector ret(vec, s );
	return ret;
}

inline Float operator*( const FloatVector& V1, const FloatVector& V2 )
{
	assert( V1.getdimension() == V2.getdimension() );
	Float ret = 0.;
	for( int p = 0; p < V1.getdimension(); p++ )
		ret += V1.getentry(p) * V2.getentry(p);
	return ret;
}	


/* Output stream notation */
inline std::ostream& operator<<( std::ostream& out, const FloatVector& V )
{
	V.print( out );
	return out;
}	
 


#endif