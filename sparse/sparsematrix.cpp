
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>
#include <cctype>

#include "sparsematrix.hpp"



SparseMatrix::SparseMatrix( int dimout, int dimin, int numentries, std::function<MatrixEntry(int)> generator )
: LinearOperator( dimout, dimin ), 
  entries(0)
{
    assert( 0 <= numentries );
    entries.reserve( numentries );
    
    for( int i = 0; i < numentries; i++ ) 
    {
      SparseMatrix::MatrixEntry entry = generator(i);
      
      assert( 0 <= entry.row && entry.row < getdimout() );
      assert( 0 <= entry.column && entry.column < getdimin() );
      
      entries.push_back( entry );
    }
    
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& entries )
: LinearOperator( dimout, dimin ), entries(entries)
{
    SparseMatrix::check();
}

SparseMatrix::SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& ent )
: LinearOperator( dimout, dimin ), entries(0)
{
    entries.reserve( ent.size() );
    std::copy( ent.begin(), ent.end(), entries.begin() );
    assert( ent.size() == entries.size() );
}

SparseMatrix::SparseMatrix( const ScalingOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getscaling() };
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( const DiagonalOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getdiagonal().at(r) };
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( const DenseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ), 
  entries( matrix.getdimout() )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        entries[ r * getdimin() + c ] = { r, c, matrix.get(r,c) };
    SparseMatrix::check();    
}





SparseMatrix::~SparseMatrix()
{
    
}

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






void SparseMatrix::reserve( int n ) const
{
    entries.reserve( n );
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
        
void SparseMatrix::setentry( int i, int r , int c, Float v )
{
//     check();
    SparseMatrix::MatrixEntry temp;
    temp.row = r;
    temp.column = c;
    temp.value = v;
    setentry( i, temp );
//     check();
}

void SparseMatrix::setentry( int i, MatrixEntry entry )
{
//     check();
    assert( 0 <= entry.row && entry.row <= getdimout() );
    assert( 0 <= entry.column && entry.column <= getdimin() );
    assert( 0 <= i && i < entries.size() );
    entries.at(i) = entry;
//     check();
}
        
        
void SparseMatrix::addentry( int r, int c, Float v )
{
//     check();
    SparseMatrix::MatrixEntry temp;
    temp.row = r;
    temp.column = c;
    temp.value = v;
    addentry( temp );
//     check();
}

void SparseMatrix::addentry( SparseMatrix::MatrixEntry entry )
{
//     check();
    assert( 0 <= entry.row && entry.row <= getdimout() );
    assert( 0 <= entry.column && entry.column <= getdimin() );
    entries.push_back( entry );
//     check();
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

    std::sort( entries.begin(), entries.end(), []( MatrixEntry a, MatrixEntry b) {
        return a.row < b.row || ( a.row == b.row && a.column < b.column );   
    });
    
//     const int N = getnumberofentries();
//     
//     for( int i = 0; i < N; i++ ) {
//     
//         MatrixEntry entry = entries.at(i);
//         int j = i;
//         while( j > 0 and ( entries.at(j-1).row > entry.row or entries.at(j-1).column > entry.column ) ) {
//             entries.at(j) = entries.at(j-1);
//             j--;
//         }
//         entries.at(j) = entry;
//     }
    
    check();
}

void SparseMatrix::sortandcompressentries() const
{
    check();
    
    sortentries();
    
    int dest = 0;
    
    for( int src = 1; src < getnumberofentries(); src++ )
    {
        assert( dest < src );
    
        if( entries[src].row == entries[dest].row and entries[src].column == entries[dest].column ) {
            entries[dest].value += entries[src].value;
        } else {
            dest++;
            entries[dest] = entries[src];
        }
        
        assert( dest <= src );
    }
    
    
    entries.resize( dest+1 );
    
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
        





SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right )
{
    left.sortandcompressentries();
    right.sortandcompressentries();
    
    int counter = 0;
    for( SparseMatrix::MatrixEntry l : left.getentries() )
    for( SparseMatrix::MatrixEntry r : left.getentries() )
        if( l.column == r.row ) 
            counter++;
            
    std::vector<SparseMatrix::MatrixEntry> new_entries;
    new_entries.reserve( counter );
    for( SparseMatrix::MatrixEntry l : left.getentries() )
    for( SparseMatrix::MatrixEntry r : left.getentries() )
        if( l.column == r.row ) 
            new_entries.push_back( { l.row, r.column, l.value * r.value } );
    
    SparseMatrix ret( left.getdimout(), right.getdimin(), new_entries );
        
    ret.sortandcompressentries();
    
    return ret;
    
}








