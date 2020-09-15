
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "densematrix.hpp"

#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/simpleoperators.hpp"
#include "../sparse/sparsematrix.hpp"


DenseMatrix::DenseMatrix( const DenseMatrix& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( new Float[ mat.getdimout() * mat.getdimin() ] ) // wierd
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = mat(r,c);
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( DenseMatrix&& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( mat.entries )
{
    assert( entries != nullptr );
    mat.entries = nullptr;
    DenseMatrix::check();
}
        
DenseMatrix& DenseMatrix::operator=( const DenseMatrix& mat )
{
    assert( entries     != nullptr );
    assert( mat.entries != nullptr );
    assert( getdimin()  == mat.getdimin()  );
    assert( getdimout() == mat.getdimout() );
    
    if( this != &mat ) {
        for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = mat(r,c);
    }
    
    DenseMatrix::check();
    return *this;
}
        
DenseMatrix& DenseMatrix::operator=( DenseMatrix&& mat ) 
{
    assert( entries     != nullptr );
    assert( mat.entries != nullptr );
    assert( getdimin()  == mat.getdimin()  );
    assert( getdimout() == mat.getdimout() );

    if( this != &mat ){
        delete[] entries;
        entries = mat.entries;
        mat.entries = nullptr;
    }

    DenseMatrix::check();
    return *this;
}
        
        
        
DenseMatrix::DenseMatrix( int dim, Float value )
: DenseMatrix( dim, dim, value )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int dim, const std::function<Float(int,int)>& generator )
: DenseMatrix( dim, dim, generator )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int dim, const std::vector<FloatVector>& coldata )
: DenseMatrix( dim, dim, coldata )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, Float value )
: LinearOperator( rows, columns ), entries( new Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = value;
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, const std::function<Float(int,int)>& generator )
: LinearOperator( rows, columns ), entries( new Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = generator(r,c);
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, const std::vector<FloatVector>& coldata )
: LinearOperator( rows, columns ), entries( new Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = coldata.at(c).at(r);
    DenseMatrix::check();
}


DenseMatrix::DenseMatrix( const ScalingOperator& scaling )
: LinearOperator( scaling.getdimout(), scaling.getdimin() ), 
  entries( new Float[ scaling.getdimout() * scaling.getdimin() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = ( ( r == c ) ? scaling.getscaling() : 0. );
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const DiagonalOperator& dia )
: LinearOperator( dia.getdimout(), dia.getdimin() ), 
  entries( new Float[ dia.getdimout() * dia.getdimin() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = ( ( r == c ) ? dia.getdiagonal().at(r) : 0. );
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const SparseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ), 
  entries( new Float[ matrix.getdimout() * matrix.getdimin() ] )
{
    assert( entries != nullptr );
    for( const SparseMatrix::MatrixEntry& entry : matrix.getentries() )
    {
        (*this)( entry.row, entry.column ) += entry.value;
    }
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const FloatVector& myvector )
: LinearOperator( myvector.getdimension(), 1 ), 
  entries( new Float[ myvector.getdimension() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < myvector.getdimension(); r++ )
    {
        (*this)( r, 0 ) = myvector[r];
    }
    DenseMatrix::check();
}
        
DenseMatrix::~DenseMatrix()
{
    // DenseMatrix::check(); // explicitly disabled, might be in moved-from state 
    
    if( entries != nullptr ){
        delete[] entries;
    }

}




void DenseMatrix::check() const 
{
    #ifdef NDEBUG
    return;
    #endif
    
    LinearOperator::check();
    assert( entries != nullptr );
}

void DenseMatrix::print( std::ostream& os ) const
{
    check();
    os << "Print Matrix " << getdimout() << "x" << getdimin() << nl;
    for( int r = 0; r < getdimout(); r++ ) {
        for( int c = 0; c < getdimin(); c++ )
            os << (*this)(r,c) << "\t";
        os << nl;
    }
    os << std::endl;
}

void DenseMatrix::printplain( std::ostream& os ) const
{
    check();
    os << getdimout() << space << getdimin() << nl;
    for( int r = 0; r < getdimout(); r++ ) {
        for( int c = 0; c < getdimin(); c++ )
            os << (*this)(r,c) << "\t";
        os << nl;
    }
    os << std::endl;
}


DenseMatrix DenseMatrix::clone() const
{
    return DenseMatrix( *this );
}




void DenseMatrix::apply( FloatVector& dest, const FloatVector& add, Float scaling ) const 
{
    add.check();
    dest.check();
    assert( dest.getdimension() == getdimout() );
    assert( add.getdimension() == getdimin() );
    
    dest.zero();
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        dest[r] += scaling * (*this)(r,c) * add[c];
    
    /*
    load r 0
    .run_outer_loop
    if r == dimout goto .end_outer_loop
      load c 0
      .run_inner_loop
      if c == dimin goto .end_inner_loop
        // perform multiplication
      inc
      goto .run_inner_loop
      .end_inner_loop
      inc r
      goto run_outer_loop
    .end_outer_loop
    */
    
    /*
    load r 0
    .run_outer_loop
        load c 0
        .run_inner_loop
        
        // perform with r and c
        
        inc c
        if c < dimin goto .run_inner_loop
    inc r
    if r < dimout goto .run_outer_loop
    */

}

Float DenseMatrix::get( int r, int c ) const
{
    return (*this)( r, c );
}

void DenseMatrix::set( int r, int c, Float v )
{
    (*this)( r, c ) = v;
}        

Float& DenseMatrix::at( int r, int c ) &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[  r * getdimin() + c ];
}

const Float& DenseMatrix::at( int r, int c ) const &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[  r * getdimin() + c ];
}

Float& DenseMatrix::operator()( int r, int c ) &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[ r * getdimin() + c ];
}

const Float& DenseMatrix::operator()( int r, int c ) const &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[ r * getdimin() + c ];
}

DenseMatrix DenseMatrix::submatrix( const IndexMap& rows, const IndexMap& columns ) const
{
    check();
    rows.check(); 
    columns.check();
    
    assert( rows.isstrictlyascending()    );
    assert( columns.isstrictlyascending() );
    
    DenseMatrix ret( rows.getSourceRange().cardinality(), columns.getSourceRange().cardinality(), 0. );
    
    for( int nr = 0; nr < ret.getdimout(); nr++ )
    for( int nc = 0; nc < ret.getdimin();  nc++ )
    {
        assert( rows.getSourceRange().contains(nr) );
        assert( columns.getSourceRange().contains(nc) );
        
        ret( nr, nc ) = (*this)( rows[nr], columns[nc] );
    }
    
    return ret;
    
//     if( getdimin() == 0 || getdimout() == 0 )
//         return *this;
//     
//     assert( rows.getSourceRange().min() == 0 );
//     assert( rows.getSourceRange().max() <= getdimout() - 1 );
//     assert( rows.getTargetRange().min() == 0 );
//     assert( rows.getTargetRange().max() == getdimout() - 1 );
//     assert( columns.getSourceRange().min() == 0 );
//     assert( columns.getSourceRange().max() <= getdimin() - 1 );
//     assert( columns.getTargetRange().min() == 0 );
//     assert( columns.getTargetRange().max() == getdimin() - 1 );
//     assert( rows.isstrictlyascending() );
//     assert( columns.isstrictlyascending() );
//     
//     DenseMatrix ret( std::max(0,rows.getSourceRange().max() - 1), 
//                     std::max(0,columns.getSourceRange().max() - 1) );
//     for( int nr = 0; nr < ret.getdimout(); nr++ )
//     for( int nc = 0; nc < ret.getdimin();  nc++ )
//         ret( nr, nc ) = (*this)( rows[nr], columns[nc] );
//     
//     ret.check();
//     return ret;
    
}




FloatVector DenseMatrix::flattencolumns() const
{
    check();
    FloatVector ret( getdimin() * getdimout() );
    for( int c = 0; c < getdimin();  c++ )
    for( int r = 0; r < getdimout(); r++ )
        ret[ c * getdimout() + r ] = get(r,c);
    return ret;
}

FloatVector DenseMatrix::flattenrows() const
{
    check();
    FloatVector ret( getdimin() * getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        ret[ r * getdimin() + c ] = get(r,c);
    return ret;
}

// std::vector<SparseMatrix::MatrixEntry> DenseMatrix::getSparseMatrixEntries( bool clean ) const
// {
//     check();
//     std::vector<SparseMatrix::MatrixEntry> ret(0);
//     ret.reserve( getdimin() * getdimout() );
//     for( int r = 0; r < getdimout(); r++ )
//     for( int c = 0; c < getdimin (); c++ )
//         if( not clean or get(r,c) != 0. )
//             ret.push_back( { r, c, get(r,c) } );
//     ret.reserve( ret.size() );
//     return ret;
// }









void DenseMatrix::zeromatrix()
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = 0.;
}

void DenseMatrix::unitmatrix()
{
    check();
    assert( getdimout() == getdimin() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        if( r == c ) 
            (*this)(r,c) = 1.;
        else
            (*this)(r,c) = 0.;
}

void DenseMatrix::randommatrix()
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = gaussrand(); // sqrt( rand() );
}

void DenseMatrix::randomintegermatrix( int min, int max )
{
    check();
    assert( min+2 <= max );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = min + rand() % (max-min);
}

void DenseMatrix::scale( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) *= s;
}

void DenseMatrix::set( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = s;
}

void DenseMatrix::add( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) += s;
}

FloatVector DenseMatrix::getrow( int r ) const 
{
    check();
    assert( 0 <= r && r < getdimout() );
    FloatVector row( getdimin() );
    for( int c = 0; c < getdimin(); c++ )
            row[c] = (*this)(r,c);
    row.check();
    return row;
}

void DenseMatrix::setrow( int r, const FloatVector& row )
{
    check();
    row.check();
    assert( 0 <= r && r < getdimout() );
    assert( row.getdimension() == getdimin() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) = row[c];
}

FloatVector DenseMatrix::getcolumn( int c ) const 
{
    check();
    assert( 0 <= c && c < getdimin() );
    FloatVector column( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        column[r] = (*this)(r,c);
    column.check();
    return column;
}

void DenseMatrix::setcolumn( int c, const FloatVector& column )
{
    check();
    column.check();
    assert( 0 <= c && c < getdimin() );
    assert( column.getdimension() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c) = column[r];
}

void DenseMatrix::addrow( int r, const FloatVector& row, Float s )
{
    check();
    row.check();
    assert( 0 <= r && r < getdimout() );
    assert( row.getdimension() == getdimin() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) += s * row[c];
}

void DenseMatrix::addcolumn( int c, const FloatVector& column, Float s )
{
    check();
    column.check();
    assert( 0 <= c && c < getdimin() );
    assert( column.getdimension() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c) += s * column[r];
}

void DenseMatrix::swaprow( int r1, int r2 )
{
    check();
    assert( 0 <= r1 && r1 < getdimout() );
    assert( 0 <= r2 && r2 < getdimout() );
    if( r1 == r2 ) return;
    for( int c = 0; c < getdimin(); c++ )
        std::swap( (*this)(r1,c), (*this)(r2,c) );
    check();
}

void DenseMatrix::swapcolumn( int c1, int c2 )
{
    check();
    assert( 0 <= c1 && c1 < getdimin() );
    assert( 0 <= c2 && c2 < getdimin() );

    if( c1 == c2 ) return;

    for( int r = 0; r < getdimout(); r++ )
        std::swap( (*this)(r,c1), (*this)(r,c2) );

    check();
}


void DenseMatrix::scalerow( int r, Float alpha )
{
    check();
    assert( 0 <= r && r < getdimout() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) *= alpha;
}

void DenseMatrix::scalecolumn( int c, Float alpha )
{
    check();
    assert( 0 <= c && c < getdimin() );
    for( int r = 0; r < getdimout(); c++ )
      (*this)(r,c) *= alpha;
}

void DenseMatrix::addrow( int r1, int r2, Float alpha )
{
    check();
    assert( 0 <= r1 && r1 < getdimout() );
    assert( 0 <= r2 && r2 < getdimout() );
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r1,c) += alpha * (*this)(r2,c);
    check();
}

void DenseMatrix::addcolumn( int c1, int c2, Float alpha )
{
    check();
    assert( 0 <= c1 && c1 < getdimin() );
    assert( 0 <= c2 && c2 < getdimin() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c1) += alpha * (*this)(r,c2);
    check();
}


        




        
void DenseMatrix::indexmapping( const IndexMap& im )
{
    check();
    im.check();
    if( im.getSourceRange().isempty() )
        return;
    assert( im.getSourceRange().min() == 0 );
    assert( im.getSourceRange().max() == getdimin() - 1 );
    assert( im.getTargetRange().min() == 0 );
    assert( im.getTargetRange().max() == getdimout() - 1 );
    zeromatrix();
    for( int c = 0; c < im.getSourceRange().max() - 1; c++ )
        (*this)( im[c], c ) = 1.;
    check();
}





        
void DenseMatrix::add( const DenseMatrix& addendum )
{
    add( 1., addendum );
}

void DenseMatrix::add( Float s, const DenseMatrix& addendum )
{
    add( 1., s, addendum );
}

void DenseMatrix::add( Float s, Float t, const DenseMatrix& addendum )
{
    check();
    addendum.check();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = s * (*this)(r,c) + t * addendum(r,c);
    check();
}
        



DenseMatrix DenseMatrix::symmetricPart() const 
{
    assert( issquare() );
    DenseMatrix ret( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret(r,c) = 0.5 * (*this)(r,c) + 0.5 * (*this)(c,r);
    return ret;
}

DenseMatrix DenseMatrix::antisymmetricPart() const 
{
    assert( issquare() );
    DenseMatrix ret( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret(r,c) = 0.5 * (*this)(r,c) - 0.5 * (*this)(c,r);
    return ret;
}






Float DenseMatrix::maxabsoluteentry() const
{
    check();
    Float ret = - std::numeric_limits<Float>::infinity();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            ret = std::max( ret, (*this)(r,c) );
    return ret;
}









bool DenseMatrix::issquare() const
{
    return getdimout() == getdimin();
}

bool DenseMatrix::issymmetric() const
{
    if( not issquare() ) return false;
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c <= r; c++ )
        if( not isabout( get(r,c), get(c,r) ) ) 
            return false;
    
    return true;
}

bool DenseMatrix::isantisymmetric() const
{
    if( not issquare() ) return false;
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c <= r; c++ )
        if( not isabout( get(r,c), -get(c,r) ) ) 
            return false;
    
    return true;
}

bool DenseMatrix::isfinite() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( not std::isfinite( get(r,c) ) ) 
            return false;
    return true;
}

bool DenseMatrix::iszero() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) != 0. ) 
            return false;
    return true;
}



bool DenseMatrix::ispositive() const 
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) <= 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnegative() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) >= 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnonnegative() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) < 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnonpositive() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) > 0. ) 
            return false;
    return true;
}

bool DenseMatrix::issmall( Float eps ) const
{
    return this->norm() < eps;
}


Float DenseMatrix::norm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += get(r,c) * get(r,c);
    return std::sqrt( ret );
}

Float DenseMatrix::maxnorm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret = std::max( ret, get(r,c) );
    return ret;
}

Float DenseMatrix::sumnorm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += absolute( get(r,c) );
    return ret;
}

Float DenseMatrix::lpnorm( Float p ) const
{
    assert( p > 1. );
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += power_numerical( absolute( get(r,c) ), p );
    return power_numerical( ret, 1. / p );
}


Float* DenseMatrix::raw()
{
    return entries;
}

const Float* DenseMatrix::raw() const
{
    return entries;
}
