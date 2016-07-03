
#include <iostream>
#include <algorithm>
#include <iterator>
#include <list>
#include <cctype>

#include "densematrix.hpp"




DenseMatrix::DenseMatrix( int rows, int columns )
: LinearOperator( rows, columns), entries( rows * columns )
{
    
}

DenseMatrix::~DenseMatrix()
{
    
}




void DenseMatrix::check() const 
{
    LinearOperator::check();
}

void DenseMatrix::print( std::ostream& os ) const
{
    os << "Print Matrix " << getdimout() << "x" << getdimin() << std::endl;
    for( int r = 0; r < getdimout(); r++ ) {
        for( int c = 0; c < getdimin(); c++ )
            os << (*this)(r,c) << "\t";
        os << std::endl;
    }
}

void DenseMatrix::applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const 
{
    assert( dest.getdimension() == getdimout() );
    assert( add.getdimension() == getdimin() );
    dest.scale( s );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        dest[r] += t * (*this)(r,c) * add[c];
}

Float DenseMatrix::get( int r, int c ) const
{
    return (*this)( r, c );
}

void DenseMatrix::set( int r, int c, Float v )
{
    (*this)( r, c ) = v;
}        

Float& DenseMatrix::operator()( int r, int c )
{
    return entries.at( r * getdimin() + c );
}

const Float& DenseMatrix::operator()( int r, int c ) const
{
    return entries.at( r * getdimin() + c );
}

DenseMatrix DenseMatrix::submatrix( const IndexMap& rows, const IndexMap& columns ) const
{
    rows.check(); columns.check();
    assert( rows.getSourceRange().getlow() == 0 );
    assert( rows.getSourceRange().gethigh() <= getdimout() - 1 );
    assert( rows.getDestRange().getlow() == 0 );
    assert( rows.getDestRange().gethigh() == getdimout() - 1 );
    assert( columns.getSourceRange().getlow() == 0 );
    assert( columns.getSourceRange().gethigh() <= getdimin() - 1 );
    assert( columns.getDestRange().getlow() == 0 );
    assert( columns.getDestRange().gethigh() == getdimin() - 1 );
    assert( rows.isstrictlyascending() );
    assert( columns.isstrictlyascending() );
    
    DenseMatrix ret( rows.getSourceRange().gethigh() - 1, columns.getSourceRange().gethigh() - 1 );
    for( int nr = 0; nr < ret.getdimout(); nr++ )
    for( int nc = 0; nc < ret.getdimin();  nc++ )
        ret( nr, nc ) = (*this)( rows[nr], columns[nc] );
    return ret;
}
        
        
void DenseMatrix::zeromatrix()
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = 0.;
}

void DenseMatrix::unitmatrix()
{
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
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = sqrt( rand() );
}

void DenseMatrix::scale( Float s )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) *= s;
}

void DenseMatrix::set( Float s )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = s;
}

FloatVector DenseMatrix::getrow( int r ) const 
{
	assert( 0 <= r && r < getdimout() );
	FloatVector row( getdimin() );
	for( int c = 0; c < getdimin(); c++ )
		row[c] = (*this)(r,c);
	return row;
}

void DenseMatrix::setrow( int r, const FloatVector& row )
{
	assert( 0 <= r && r < getdimout() );
	assert( row.getdimension() == getdimin() );
	for( int c = 0; c < getdimin(); c++ )
		(*this)(r,c) = row[c];
}

FloatVector DenseMatrix::getcolumn( int c ) const 
{
	assert( 0 <= c && c < getdimin() );
	FloatVector column( getdimout() );
	for( int r = 0; r < getdimout(); r++ )
		column[r] = (*this)(r,c);
	return column;
}

void DenseMatrix::setcolumn( int c, const FloatVector& column )
{
	assert( 0 <= c && c < getdimin() );
	assert( column.getdimension() == getdimout() );
	for( int r = 0; r < getdimout(); r++ )
		(*this)(r,c) = column[r];
}

void DenseMatrix::indexmapping( const IndexMap& im )
{
    assert( im.getSourceRange().getlow() == 0 );
    assert( im.getSourceRange().gethigh() <= getdimin() - 1 );
    assert( im.getDestRange().getlow() == 0 );
    assert( im.getDestRange().gethigh() == getdimout() - 1 );
    zeromatrix();
    for( int c = 0; c < im.getSourceRange().gethigh() - 1; c++ )
        (*this)( c, im[c] ) = 1.;
}

DenseMatrix DenseMatrix::transpose() const
{
    DenseMatrix ret( getdimin(), getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        ret.entries.at( c * getdimout() + r ) = entries.at( r * getdimin() + c );
    return ret;
}


// DenseMatrix DenseMatrix::inverse() const
// {
    // int d = getdimout();
    
    // DenseMatrix l(d,d), u(d,d), li(d,d), ui(d,d), ret(d,d);
    // Float det;
    
    // gaussfactorization( l, u, li, ui, det );
    
    // ret = ui * li;
    
    // return ret;
// }
        

// void DenseMatrix::gaussfactorization( 
        // DenseMatrix& Lower, DenseMatrix& Upper,
        // DenseMatrix& LowerInv, DenseMatrix& UpperInv,
        // Float& det ) const
// {

    // assert( this->issquare() && Lower.issquare() && Upper.issquare() );
    // assert( this->getdimout() == Lower.getdimout() && Lower.getdimout() == Upper.getdimout() );
    // Lower.set( 47.11 );
    // Upper.set( 47.11 );
    
    // const int n = this->getdimout();
    
    // /* Construct lower and upper matrix */
    
    // for( int i = 0; i < n; i++ )
    // {
        
        // for( int j = 0; j < n; j++ )
            // if (j < i)
                // Lower( j, i ) = 0.;
            // else
            // {
                // Lower( j, i ) = (*this)( j, i );
                // for( int k = 0; k < i; k++)
                    // Lower( j, i ) = Lower( j, i ) - Lower( j, k ) * Upper( k, i );
            // }

        // for( int j = 0; j < n; j++ )
            // if (j < i)
                // Upper( i, j ) = 0.;
            // else if ( j == i )
                // Upper( i, j ) = 1.;
            // else
            // {
                // Upper( i, j ) = (*this)( i, j ) / Lower( i, i );
                // for( int k = 0; k < i; k++ )
                    // Upper( i, j ) = Upper( i, j ) - ( ( Lower( i, k ) * Upper( k, j ) ) / Lower( i, i ));
            // }

    // }
    
    // /* Determinant */
    
    // det = 1.;
    // for( int i = 0; i < n; i++ )
        // det *= Upper(i,i) * Lower(i,i);
    
    // /* Construct inverses of the triangular matrices */
    
    // assert( LowerInv.issquare() && UpperInv.issquare() );
    // assert( this->getdimout() == LowerInv.getdimout() && LowerInv.getdimout() == UpperInv.getdimout() );
    // LowerInv.set( 47.11 );
    // UpperInv.set( 47.11 );
    
    // // TODO: Compute the Inverses via Back Substitution

// }
        
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
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = s * (*this)(r,c) + t * addendum(r,c);
}
        
DenseMatrix DenseMatrix::MatrixMult( const DenseMatrix& left, const DenseMatrix& right )
{
    int lin = left.getdimin();
    int lout = left.getdimout();
    int rin = right.getdimin();
    int rout = right.getdimout();
    
    assert( lin == rout );
    
    DenseMatrix ret( lout, rin );
    ret.zeromatrix();
    
    for( int lo = 0; lo < lout; lo++ )
        for( int ri = 0; ri < rin; ri++ )
            for( int m = 0; m < rout; m++ )
                ret( lo, ri ) += left( lo, m ) * right( m, ri );
    
    return ret;
    
}

bool DenseMatrix::issquare() const
{
    return getdimout() == getdimin();
}