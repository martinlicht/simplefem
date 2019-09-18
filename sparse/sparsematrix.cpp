
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
    #ifdef NDEBUG
    return;
    #endif
    
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
		

		
bool SparseMatrix::is_sorted( SparseMatrix::MatrixEntrySorting manner ) const
{
    for( int i = 1; i < entries.size(); i++ ) {

        const auto& a = entries[i-1];
        const auto& b = entries[  i];

        if( manner == MatrixEntrySorting::rowwise ){
            if( a.row    > b.row    or ( a.row == b.row and a.column > b.column ) )
                return false;
        }else{
            if( a.column > b.column or ( a.column == b.column and a.row > b.row ) )
                return false;
        }

    }

    return true;
    
}

void SparseMatrix::sortentries( SparseMatrix::MatrixEntrySorting manner ) const
{
    check();

    if( manner == MatrixEntrySorting::rowwise )
        std::sort( entries.begin(), entries.end(), []( MatrixEntry a, MatrixEntry b) {
            return a.row < b.row || ( a.row == b.row && a.column < b.column );   
        });
    else 
        std::sort( entries.begin(), entries.end(), []( MatrixEntry a, MatrixEntry b) {
            return a.column < b.column || ( a.column == b.column && a.row < b.row );   
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
    
    assert( is_sorted( manner ) );

    check();
}

void SparseMatrix::sortandcompressentries( SparseMatrix::MatrixEntrySorting manner ) const
{
    check();
    
    sortentries( manner );
    
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
    
    auto newentries = entries;

    for( auto& newentry : newentries )
        std::swap( newentry.row, newentry.column );

    return SparseMatrix( getdimin(), getdimout(), newentries );

    // SparseMatrix ret = SparseMatrix( getdimout(), getdimin() );
    
    // for( const MatrixEntry& me : entries ) {
    //     ret.entries.push_back( { me.column, me.row, me.value } );
    // }
    
    // return ret;
}
        







SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right )
{

    std::cout << "--- SparseMatrix Product" << std::endl;
    std::cout << "--- Sort and compress" << std::endl;
    
    left.sortandcompressentries( SparseMatrix::MatrixEntrySorting::columnwise );
    right.sortandcompressentries( SparseMatrix::MatrixEntrySorting::rowwise );

    std::cout << "--- Counting" << std::endl;
    
    int counter = 0;
    
    {

        int lbase = 0; int rbase = 0;
        while( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() ){
            assert( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() );
            while( lbase < left.getnumberofentries() and left.getentry(lbase).column < right.getentry(rbase).row ) lbase++;
            if( lbase == left.getnumberofentries() ) break;
            while( rbase < right.getnumberofentries() and right.getentry(rbase).row < left.getentry(lbase).column ) rbase++;
            if( rbase == right.getnumberofentries() ) break;
            assert( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() );
            assert( left.getentry(lbase).column == right.getentry(rbase).row );
            int index = left.getentry(lbase).column;
            int lnext = lbase; int rnext = rbase;
            while( lnext < left.getnumberofentries() and left.getentry(lnext).column == index ) lnext++;
            while( rnext < right.getnumberofentries() and right.getentry(rnext).row == index  ) rnext++;
            counter = counter + ( lnext - lbase ) * ( rnext - rbase );
            lbase = lnext; rbase = rnext;
        }

    }

    std::cout << "--- Assemble" << std::endl;
    
    std::vector<SparseMatrix::MatrixEntry> new_entries;
    new_entries.reserve( counter );
    
    {
        
        int lbase = 0; int rbase = 0;
        while( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() ){
            assert( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() );
            while( lbase < left.getnumberofentries() and left.getentry(lbase).column < right.getentry(rbase).row ) lbase++;
            if( lbase == left.getnumberofentries() ) break;
            while( rbase < right.getnumberofentries() and right.getentry(rbase).row < left.getentry(lbase).column ) rbase++;
            if( rbase == right.getnumberofentries() ) break;
            assert( lbase < left.getnumberofentries() and rbase < right.getnumberofentries() );
            assert( left.getentry(lbase).column == right.getentry(rbase).row );
            int index = left.getentry(lbase).column;
            int lnext = lbase; int rnext = rbase;
            while( lnext < left.getnumberofentries() and left.getentry(lnext).column == index ) lnext++;
            while( rnext < right.getnumberofentries() and right.getentry(rnext).row == index  ) rnext++;
            for( int lcurr = lbase; lcurr < lnext; lcurr++ )
            for( int rcurr = rbase; rcurr < rnext; rcurr++ )
                new_entries.push_back( { 
                                        left.getentry(lcurr).row, 
                                        right.getentry(rcurr).column, 
                                        left.getentry(lcurr).value * right.getentry(rcurr).value
                                       } );
            lbase = lnext; rbase = rnext;
        }

    }

    assert( new_entries.size() == counter );

    std::cout << "--- Construct" << std::endl;
    SparseMatrix ret( left.getdimout(), right.getdimin(), new_entries );
        
    // std::cout << "--- Sort and compress again" << std::endl;
    // ret.sortandcompressentries();
    
    return ret;
    
}







SparseMatrix SparseMatrixMultiplication( const SparseMatrix& left, const SparseMatrix& right )
{

    std::cout << "--- SparseMatrix Product" << std::endl;
    std::cout << "--- Sort and compress" << std::endl;
    
    left.sortandcompressentries( SparseMatrix::MatrixEntrySorting::columnwise );
    right.sortandcompressentries( SparseMatrix::MatrixEntrySorting::rowwise );

    std::cout << "--- Counting" << std::endl;
    
    int counter = 0;
    for( SparseMatrix::MatrixEntry l : left.getentries()  )
    for( SparseMatrix::MatrixEntry r : right.getentries() )
        if( l.column == r.row ) 
            counter++;

    std::cout << "--- Assemble" << std::endl;
    
    std::vector<SparseMatrix::MatrixEntry> new_entries;
    new_entries.reserve( counter );
    for( SparseMatrix::MatrixEntry l : left.getentries()  )
    for( SparseMatrix::MatrixEntry r : right.getentries() )
        if( l.column == r.row ) 
            new_entries.push_back( { l.row, r.column, l.value * r.value } );
    
    assert( new_entries.size() == counter );

    std::cout << "--- Construct" << std::endl;
    SparseMatrix ret( left.getdimout(), right.getdimin(), new_entries );
        
    // std::cout << "--- Sort and compress again" << std::endl;
    // ret.sortandcompressentries();
    
    return ret;
    
}








