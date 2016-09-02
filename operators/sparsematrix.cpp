
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "sparsematrix.hpp"


SparseMatrix::SparseMatrix( int dimout, int dimin )
: LinearOperator(dimout,dimin)
{
    check();
}

SparseMatrix::SparseMatrix( const ScalingOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() )
{
    for( int r = 0; r < getdimout(); r++ )
        entries.push_back( { r, r, matrix.getscaling() } );
    check();    
}

SparseMatrix::SparseMatrix( const DiagonalOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() )
{
    for( int r = 0; r < getdimout(); r++ )
        entries.push_back( { r, r, matrix.getdiagonal().at(r) } );
    check();    
}

SparseMatrix::SparseMatrix( const DenseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        entries.push_back( { r, c, matrix.get(r,c) } );
    check();    
}



SparseMatrix::~SparseMatrix()
{}
	
void SparseMatrix::check() const
{
	LinearOperator::check();
    for( const MatrixEntry& me : entries ) {
        assert( 0 <= me.row && me.row < getdimout() );
        assert( 0 <= me.column && me.column < getdimin() );
    }
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
	SparseMatrix::MatrixEntry temp;
	temp.row = r;
	temp.column = c;
	temp.value = v;
	addentry( temp );
}

void SparseMatrix::addentry( SparseMatrix::MatrixEntry entry )
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

const std::list<SparseMatrix::MatrixEntry>& SparseMatrix::getentries() const
// const std::list<MatrixEntry>& SparseMatrix::getentries() const
{
	return entries;
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
	
void SparseMatrix::compressentries() const
{
	sortentries();
	
	auto mergeMatrixEntry = []( const SparseMatrix::MatrixEntry& left, 
					const SparseMatrix::MatrixEntry& right)
	-> SparseMatrix::MatrixEntry {
	 SparseMatrix::MatrixEntry ret;
	 assert( left.row == right.row && left.column == right.column );
	 ret.row = left.row; ret.column = left.column;
	 ret.value = left.value + right.value;
     return ret; 
	};
						
	mergeelementsinsortedlist<SparseMatrix::MatrixEntry>
	( const_cast<std::list<MatrixEntry>&>(entries), 
	mergeMatrixEntry, compareMatrixEntry );
	
}
	
