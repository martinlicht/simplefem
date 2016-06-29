#ifndef INCLUDEGUARD_SPARSEMATRIX
#define INCLUDEGUARD_SPARSEMATRIX

#include <list>
#include "../basic.hpp"
#include "linearoperator.hpp"


	
class SparseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

	public:
	
	struct MatrixEntry
	{
		int row;
		int column;
		Float value;
	};
	
	static bool compareMatrixEntry(const MatrixEntry& first, const MatrixEntry& second) 
	{
            if( first.row < second.row )
                return true;
            else if( first.row == second.row && first.column < second.column )
                return true;
            else 
                return false;
	}

	explicit SparseMatrix(int,int);
	virtual ~SparseMatrix();
	
	virtual void check() const override;
	virtual void print( std::ostream& ) const override;
	
	void addentry( int, int, Float );
    void addentry( MatrixEntry );
	void clearentries();
	
	int getnumberofentries() const;
	void sortentries() const;
	
	virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;
  
	
	private:
	
	std::list<MatrixEntry> entries;
	
};
  
  

#endif