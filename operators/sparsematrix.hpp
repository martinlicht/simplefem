#ifndef INCLUDEGUARD_SPARSEMATRIX
#define INCLUDEGUARD_SPARSEMATRIX

#include <list>
#include "../basic.hpp"
#include "linearoperator.hpp"


struct MatrixEntry
	{
		int row;
		int column;
		Float value;
	};

	
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
		return first.row > second.row || first.column > second.column;
	}

	SparseMatrix(int,int);
	virtual ~SparseMatrix();
	
	virtual void check() const;
	virtual void print( std::ostream& ) const;
	
	void addentry( int, int, Float );
    void addentry( MatrixEntry );
	void clearentries();
	
	int getnumberofentries() const;
	void sortentries() const;
	
	virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const;
  
	
	private:
	
	std::list<MatrixEntry> entries;
	
};
  
  

#endif