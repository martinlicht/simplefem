
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "sparsematrix.hpp"


SparseMatrix::SparseMatrix( int dimout, int dimin )
: LinearOperator(dimout,dimin)
{}


SparseMatrix::~SparseMatrix()
{}
	
void SparseMatrix::check() const
{
	LinearOperator::check();
}

void SparseMatrix::print( std::ostream& os ) const
{
	os << "SparseMatrix Entries:" << std::endl;
	for( const MatrixEntry& entry : entries )
		os << entry.row << " " << entry.column << " : " << entry.value << std::endl; 
}
	

void SparseMatrix::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
	assert( getdimout() == dest.getdimension() );
	assert( getdimin() == add.getdimension() );
	
	dest.scale( s );
        
        for( const MatrixEntry& rcv : entries )
		dest.setentry( rcv.row, dest.getentry( rcv.row ) + t * rcv.value * add.getentry( rcv.column ) );
}


void SparseMatrix::addentry( int r, int c, Float v )
{
	MatrixEntry temp;
	temp.row = r;
	temp.column = c;
	temp.value = v;
	addentry( temp );
}

void SparseMatrix::addentry( MatrixEntry entry )
{
	assert( 0 <= entry.row && entry.row <= getdimout() );
	assert( 0 <= entry.column && entry.column <= getdimin() );
	entries.push_back( entry );
}

int SparseMatrix::getnumberofentries() const 
{
	return entries.size();
}
	
void SparseMatrix::clearentries()
{
	entries.clear();
}
		

		
void SparseMatrix::sortentries() const
{
	
//         auto& casted_entries = const_cast<std::list<MatrixEntry>&>(entries);
        
//         std::list<MatrixEntry> test = entries;
//         
//         auto comp = []( const SparseMatrix::MatrixEntry& left, const SparseMatrix::MatrixEntry& right) -> bool
//         { 
//                         return left.row < right.row || left.column < right.column || left.value < right.value;
//         };
        
	const_cast<std::list<MatrixEntry>&>(entries).sort( compareMatrixEntry );
	
}
	
	
