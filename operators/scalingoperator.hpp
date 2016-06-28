#ifndef INCLUDEGUARD_SCALINGOPERATOR
#define INCLUDEGUARD_SCALINGOPERATOR

#include "../basic.hpp"
#include "linearoperator.hpp"
// #include "sparsematrix.hpp"


class ScalingOperator:
public LinearOperator /* every matrix is a linear operator */
{

	public:
		
		ScalingOperator( int, Float s );
		virtual ~ScalingOperator();
		
		virtual void check() const override	;
		// SparseMatrix cloneSparseMatrix() const;
		
		Float getscaling();
		void setscaling( Float s );
		
		virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const;
		
	private:
	
		Float scaling;
};
  
  

#endif