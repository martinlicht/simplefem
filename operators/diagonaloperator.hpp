#ifndef INCLUDEGUARD_DIAGONALOPERATOR
#define INCLUDEGUARD_DIAGONALOPERATOR

#include "../basic.hpp"
#include "linearoperator.hpp"



class DiagonalOperator:
public LinearOperator /* every matrix is a linear operator */
{

	public:
		
		explicit DiagonalOperator( int, const FloatVector& dia );
		virtual ~DiagonalOperator();
		
		virtual void check() const override;
		virtual void print( std::ostream& ) const override;
		
		FloatVector& getdiagonal();
		const FloatVector& getdiagonal() const;
		
		virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;
		
	private:
	
		FloatVector diagonal;
		
};
  
  

#endif