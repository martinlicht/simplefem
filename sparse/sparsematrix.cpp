
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>
#include <cctype>

#include "sparsematrix.hpp"


SparseMatrix::SparseMatrix( int dimout, int dimin )
: LinearOperator(dimout,dimin)
{
    check();
}

SparseMatrix::SparseMatrix( int dimout, int dimin, int numentries )
: LinearOperator(dimout,dimin), 
  entries(numentries)
{
    check();
}

SparseMatrix::SparseMatrix( const ScalingOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getscaling() };
    check();    
}

SparseMatrix::SparseMatrix( const DiagonalOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getdiagonal().at(r) };
    check();    
}

SparseMatrix::SparseMatrix( const DenseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ), 
  entries( matrix.getdimout() )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        entries[ r * getdimin() + c ] = { r, c, matrix.get(r,c) };
    check();    
}

SparseMatrix::SparseMatrix( int dimout, int dimin, int numentries, std::function<MatrixEntry(int)> generator )
: LinearOperator( dimout, dimin ), 
  entries(0)
{
    assert( 0 <= numentries );
    entries.reserve( numentries );
    
    for( int i = 0; i < numentries; i++ ) 
    {
      entries.push_back( generator(i) );
    }
    
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
    os << std::endl;
}

void SparseMatrix::printplain( std::ostream& os ) const
{
    for( const MatrixEntry& entry : entries )
        os << entry.row << " " << entry.column << " " << entry.value << nl;
    os << std::endl;
}


FloatVector SparseMatrix::apply( const FloatVector& add, Float scaling ) const 
{
    check();
    add.check();
    assert( getdimin() == add.getdimension() );
    
    FloatVector ret( getdimout(), 0. );
    for( const MatrixEntry& rcv : entries )
        ret.setentry( rcv.row, ret.getentry( rcv.row ) + scaling * rcv.value * add.getentry( rcv.column ) );

    return ret;
}



const SparseMatrix::MatrixEntry& SparseMatrix::getentry( int i ) const
{
    assert( 0 <= i && i < getnumberofentries() );
    return entries[i];
}

SparseMatrix::MatrixEntry& SparseMatrix::getentry( int i )
{
    assert( 0 <= i && i < getnumberofentries() );
    return entries[i];
}
        
void SparseMatrix::addentry( int r, int c, Float v )
{
    check();
    SparseMatrix::MatrixEntry temp;
    temp.row = r;
    temp.column = c;
    temp.value = v;
    addentry( temp );
    check();
}

void SparseMatrix::addentry( SparseMatrix::MatrixEntry entry )
{
    check();
    assert( 0 <= entry.row && entry.row <= getdimout() );
    assert( 0 <= entry.column && entry.column <= getdimin() );
    entries.push_back( entry );
    check();
}

int SparseMatrix::getnumberofentries() const 
{
    return entries.size();
}
	
void SparseMatrix::clearentries()
{
    check();
    entries.clear();
    check();
}

const std::vector<SparseMatrix::MatrixEntry>& SparseMatrix::getentries() const
{
    return entries;
}
		

		
void SparseMatrix::sortentries() const
{
   check();
//         auto& casted_entries = const_cast<std::list<MatrixEntry>&>(entries);
        
//         std::list<MatrixEntry> test = entries;
//         
//         auto comp = []( const SparseMatrix::MatrixEntry& left, const SparseMatrix::MatrixEntry& right) -> bool
//         { 
//                         return left.row < right.row || left.column < right.column || left.value < right.value;
//         };
        
//     const_cast<std::vector<MatrixEntry>&>(entries).sort( compareMatrixEntry );
   
   // TODO: Use STL for sorting the entries of a std::vector
    
    check();
}

void SparseMatrix::compressentries() const
{
    check();
    
    sortentries();
    
    // TODO: 
    // Develop compression of sparse matrix
    // when the entries are in a std::vector 
    
//     auto mergeMatrixEntry 
//     =
//     []( const SparseMatrix::MatrixEntry& left, 
//         const SparseMatrix::MatrixEntry& right )
//       -> SparseMatrix::MatrixEntry
//       {
//         SparseMatrix::MatrixEntry ret;
//         assert( left.row == right.row && left.column == right.column );
//         ret.row = left.row; ret.column = left.column;
//         ret.value = left.value + right.value;
//         return ret; 
//       };
//                                             
//     mergeelementsinsortedlist<SparseMatrix::MatrixEntry>
//     ( const_cast<std::list<MatrixEntry>&>(entries), 
//     mergeMatrixEntry, compareMatrixEntry );
    
    check();  
}





SparseMatrix SparseMatrix::getTranspose() const 
{
    
    SparseMatrix ret = SparseMatrix( getdimout(), getdimin() );
    
    for( const MatrixEntry& me : entries ) {
        ret.entries.push_back( { me.column, me.row, me.value } );
    }
    
    return ret;
}
        













