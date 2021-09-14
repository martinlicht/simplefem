

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Conjugate Residual Method" );

int main()
{
	LOG << "Unit Test: " << TestName << endl;
	
	{
		
		LOG << "First Something Simple with the CRM" << endl;
		
		ScalingOperator A( 10, 3.141 );
		ConjugateResidualMethod CRM(A);

		FloatVector rhs(10), x(10);
		x.random(); rhs.zero();
		
		LOG << x << endl;
		CRM.solve( x, rhs );
		// LOG << x << endl;
	
	}
	

	{
		
		LOG << "First Something Simple with the Precon-CRM" << endl;
		
		ScalingOperator A( 10, Constants::pi );
		ScalingOperator M( 10, Constants::euler );
		PreconditionedConjugateResidualMethod CRM(A,M);

		FloatVector rhs(10), x(10);
		x.random(); rhs.zero();
		
		LOG << x << endl;
		CRM.solve( x, rhs );
		// LOG << x << endl;
	
	}
	

	if(true){
		
		LOG << "Now something more complicated." << endl;
		
		int dimension = 100;
        
		FloatVector x( dimension );
		for( int p = 0; p < dimension; p++ )
			x.setentry( p, 3. + p * 5. );
		
		SparseMatrix A( dimension, dimension );
		A.reserve( 3 * dimension );
        
        for( int i = 0; i < dimension; i++ ){
			if( i-1 >= 0 ) 
				A.addentry( i, i-1, 1.25 );
			if( i+1 < dimension ) 
				A.addentry( i, i+1, 1.25 );
			A.addentry( i, i, 2.51 );
		}
		// A.sortentries();
		
		SparseMatrix M( dimension, dimension );
		M.reserve( dimension );
        
        for( int i = 0; i < dimension; i++ ){
			M.addentry( i, i, 1./2.51 );
		}
		// A.sortentries();
		
        LOG << "Compute stuff." << endl;
		
		FloatVector b = A * x;
		
        {
            FloatVector y( dimension );
            srand(0);
			y.random();
            
            ConjugateResidualMethod CRM(A);
            CRM.max_iteration_count = dimension;
            CRM.threshold = 1e-20;
            
            timestamp start, end;
            start = gettimestamp();
            CRM.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
        {
            FloatVector y( dimension );
            srand(0);
			y.random();
            
            PreconditionedConjugateResidualMethod PCRM(A,M);
            PCRM.max_iteration_count = dimension;
            PCRM.threshold = 1e-20;
            
            timestamp start, end;
            start = gettimestamp();
            PCRM.solve(y,b);
            end = gettimestamp();
            LOG << timestamp2measurement( end - start ) << endl;
        }
        
		
	}
	
	LOG << "Finished Unit Test: " << TestName << endl;

	return 0;
}
