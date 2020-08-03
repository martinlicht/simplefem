

#include <cmath>
#include <algorithm>
#include <limits>


#include "../basic.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../dense/cholesky.hpp"
#include "mesh.hpp"


#ifdef NDEBUG
const int Mesh::nullindex = std::numeric_limits<int>::max(); 
#else
const int Mesh::nullindex = -17; 
#endif


Mesh::Mesh( int inner, int outer )
: innerdimension(inner), 
  outerdimension(outer),
  coordinates(outer,0)
{
  
  /* Build up the static auxiliary data */
  auxdata.resize( innerdimension+1 );
  for( int sup = 0; sup <= innerdimension; sup++ ) 
  {
  
    auxdata[sup].resize( sup+1 );
    
    for( int sub = 0; sub <= sup; sub++ )
    {
        IndexRange from( 0, sub );
        IndexRange to( 0, sup );
        std::vector<IndexMap> sigmas = generateSigmas( from, to );
        auxdata[ sup ][ sub ] = sigmas;
//         auxdata[ std::pair<int,int>(sup,sub) ] = sigmas;
    }
  
  }

  Mesh::check();
  
}

Mesh::~Mesh()
{
  Mesh::check();
}


int Mesh::getinnerdimension() const
{
  return innerdimension;
}

int Mesh::getouterdimension() const
{
  return outerdimension;
}

Coordinates& Mesh::getcoordinates()
{
  return coordinates;
}

const Coordinates& Mesh::getcoordinates() const
{
  return coordinates;
}


void Mesh::check() const 
{
  
  #ifdef NDEBUG
  return;
  #endif
    
    // check dimension 
  assert( innerdimension >= 0 );
  assert( outerdimension >= 0 );
  assert( outerdimension == coordinates.getdimension() );
  
  // the vertices and volumes must be counted 
  // assert( dimension_counted( innerdimension ) );
  // assert( dimension_counted( 0 ) );
  
  // the counting of the vertices must agree 
  // assert( count_simplices(0) == coordinates.getnumber() );
  
  /* * Data integrity
     * 
     * - for each simplex, check that the subsimplices are non-null and unique
     * - for each supersimplex, check that the supersimplices contain that simplex
     * - morphism property of inclusion 
     *   
     */
  
  // TODO: Implement 
  
  // TODO: Develop indexing rules 
  
}





int Mesh::index_from_pair( int sup, int sub ) const
{
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  return sub + (sup+1) * sup / 2;
}

void Mesh::index_to_pair( int index, int& sup, int& sub ) const
{
  // FIXME: improve this incredibly ineffecient computation.
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  for( int _sup = 1; _sup <= innerdimension; _sup++ )
  for( int _sub = 0; _sub <=            sup; _sub++ )
    if( index == index_from_pair( _sup, _sub ) )
    {
      sup = _sup; sub = _sub; return;
    }
}

int Mesh::count_subsimplices( int sup, int sub ) const 
{
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  return binomial_integer( sup + 1, sub + 1);
}




 /*
  * 
  * Accessing subsimplices 
  *
  */

bool Mesh::is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im.rangecontains( cellsub );
}

int  Mesh::get_subsimplex_index( int sup, int sub, int cellsup, int cellsub ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im.preimageof( cellsub );
}

int Mesh::get_subsimplex( int sup, int sub, int cellsup, int localindex ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im[ localindex ];  
}





 /*
  * 
  * Accessing supersimplices 
  *
  */

bool Mesh::is_supersimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  
  if( subsimplices_listed( sup, sub ) ) {
    
    return is_subsimplex( sup, sub, cellsup, cellsub );
    
  } else {
    
    assert( supersimplices_listed( sup, sub ) );
    std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
    return std::find( parents.begin(), parents.end(), cellsup ) != parents.end();
    
  }
  
}

int Mesh::get_firstparent_of_subsimplex( int sup, int sub, int cellsub ) const
{
  assert( supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  return parents[0];
}

int Mesh::get_nextparent_of_subsimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  assert( supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  it++;
  if( it == parents.end() )
    return nullindex;
  else
    return *it;
}

int Mesh::get_nextparent_by_localindex( int sup, int sub, int cellsup, int localindex ) const
{
  assert( supersimplices_listed( sup, sub ) );
  assert( subsimplices_listed( sup, sub ) );
  int cellsub = getsubsimplices( sup, sub, cellsup )[localindex];
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  it++;
  if( it == parents.end() )
    return nullindex;
  else
    return *it;
}

int Mesh::get_index_of_supersimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  assert( supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  return it - parents.begin();
}

int Mesh::get_supersimplex_by_index( int sup, int sub, int cellsub, int parentindex ) const
{
  assert( supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  return parents[ parentindex ];
}













void Mesh::set_flags( int dim, SimplexFlag flag )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( dimension_counted( dim ) );
    for( int s = 0; s < count_simplices(dim); s++ )
        set_flag( dim, s, flag );
}

const std::vector<SimplexFlag> Mesh::get_flags( int dim ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( dimension_counted( dim ) );
    std::vector<SimplexFlag> flags( count_simplices(dim), SimplexFlagInvalid );
    for( int s = 0; s < count_simplices(dim); s++ )
        flags[s] = get_flag( dim, s );
    return flags;
}

void Mesh::set_flags( int dim, std::vector<SimplexFlag> flags )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( dimension_counted( dim ) );
    assert( flags.size() == count_simplices( dim ) );
    for( int s = 0; s < count_simplices(dim); s++ )
        set_flag( dim, s, flags[s] );
}

void Mesh::automatic_dirichlet_flags()
{
    const int full = getinnerdimension();
    
    assert( dimension_counted(full-1) );
    assert( supersimplices_listed(full,full-1) );
    
    for( int d = 0; d <= getinnerdimension(); d++ )
        set_flags( d, SimplexFlagNull );
    
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_firstparent_of_subsimplex( full, full-1, s ) == nullindex || get_nextparent_of_subsimplex( full, full-1, get_firstparent_of_subsimplex( full, full-1, s ), s ) == nullindex )
            set_flag( full-1, s, SimplexFlagDirichlet );
        
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_flag( full-1, s ) == SimplexFlagDirichlet )
            for( int d = 0; d < full-1; d++ )
                for( int subindex = 0; subindex < count_subsimplices( full-1, d ); subindex++ ) 
                    set_flag( d, get_subsimplex( full-1, d, s, subindex ), SimplexFlagDirichlet );

}


void Mesh::check_dirichlet_flags()
{
    const int full = getinnerdimension();
    
    assert( dimension_counted(full-1) );
    assert( supersimplices_listed(full,full-1) );
    
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_firstparent_of_subsimplex( full, full-1, s ) == nullindex || get_nextparent_of_subsimplex( full, full-1, get_firstparent_of_subsimplex( full, full-1, s ), s ) == nullindex )
            assert( get_flag( full-1, s ) == SimplexFlagDirichlet );
        else 
            assert( get_flag( full-1, s ) == SimplexFlagNull      );
    
    for( int d = 0; d < full-1; d++ )
        for( int sub = 0; sub < count_simplices(d); sub++ )
            if( get_flag( d, sub ) == SimplexFlagDirichlet ) {
                bool found = false;
                for( int sup = get_firstparent_of_subsimplex( full-1, d, sub ); sup != nullindex; sup = get_nextparent_of_subsimplex( full-1, d, sup, sub ) )
                    found = found || ( get_flag(full-1,sup) == SimplexFlagDirichlet );
                assert( found );
            }
            
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_flag( full-1, s ) == SimplexFlagDirichlet )
            for( int d = 0; d < full-1; d++ )
                for( int subindex = 0; subindex < count_subsimplices( full-1, d ); subindex++ ) 
                    assert( get_flag( d, get_subsimplex( full-1, d, s, subindex ) ) == SimplexFlagDirichlet );
    
    
}















Float Mesh::getDiameter( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix vcm = getVertexCoordinateMatrix( dim, index );
    
    DenseMatrix dist(dim+1,dim+1);
    for( int i = 0; i <= dim; i++ )
    for( int j = 0; j <= dim; j++ )
        dist(i,j) = ( vcm.getcolumn(i) - vcm.getcolumn(j) ).norm();
    
//     assert( dist.isnonnegative() && dist.isfinite() );
    
    return dist.maxabsoluteentry();
}

Float Mesh::getMeasure( int dim, int index ) const 
{
    assert( 0 <= dim   && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix Jac = getTransformationJacobian( dim, index );
    
    DenseMatrix temp = Transpose( Jac ) * Jac;
    
    return std::sqrt(absolute(Determinant(temp))) / factorial_numerical( getinnerdimension() );
}

Float Mesh::getShapemeasure( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    return power_numerical( getDiameter( dim, index ), dim ) / getMeasure( dim, index );
}

Float Mesh::getShapemeasure( int dim ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    
    Float shapemeas = 0.;
    for( int s = 0; s < count_simplices(dim); s++ )
        shapemeas = std::max( shapemeas, power_numerical( getDiameter( dim, s ), dim ) / getMeasure( dim, s ) );
    return shapemeas;
}

Float Mesh::getShapemeasure() const 
{
    return getShapemeasure( getinnerdimension() ); 
}




DenseMatrix Mesh::getVertexCoordinateMatrix( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix ret( getouterdimension(), dim+1 );
    
    for( int v = 0; v <= dim; v++ )
    for( int d = 0; d < getouterdimension(); d++ )
        ret( d, v ) = coordinates.getdata( get_subsimplex( dim, 0, index, v ), d );
    
    return ret;
}


DenseMatrix Mesh::getTransformationJacobian( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix ret( getouterdimension(), dim );
    
    DenseMatrix vcm = getVertexCoordinateMatrix( dim, index );
    
    for( int v = 0; v < dim; v++ )
    for( int c = 0; c < getouterdimension(); c++ )
        ret( c, v ) = vcm( c, v+1 ) - vcm( c, 0 );
    
    return ret;
}


DenseMatrix Mesh::getGradientProductMatrix( int dim, int index ) const // TODO: is this the vector mass matrix?
{
    assert( 0 <= dim   && dim   <= getinnerdimension() );
    assert( 0 <= index && index <  count_simplices(dim) );
    
    DenseMatrix multiplier( dim, dim+1, 0. );
    for( int i = 0; i <  dim; i++ ) {
        multiplier(i,i+1) =  1.;
        multiplier(i,  0) = -1.;
    }
    
    // D^-1 D^-t = ( D^t D )^-1
    DenseMatrix Jac    = getTransformationJacobian( dim, index );
    DenseMatrix middle = Inverse( Transpose(Jac) * Jac );// Transpose(Jac) * Jac;//TODO: Understand
    
    return Transpose(multiplier) * middle * multiplier;
}
        
DenseMatrix Mesh::getGradientProductMatrixRightFactor( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix multiplier( dim, dim+1, 0. );
    for( int i = 0; i <  dim; i++ ) {
        multiplier(i,i+1) =  1.;
        multiplier(i,  0) = -1.;
    }
    
    // D^-1 D^-t = ( D^t D )^-1
    DenseMatrix Jac    = getTransformationJacobian( dim, index );
    DenseMatrix middle = Inverse( Transpose(Jac) * Jac ); //Transpose(Jac) * Jac; //TODO: Understand
    
    DenseMatrix middle_rightfactor = Transpose( CholeskyDecomposition( middle ) ); 

    assert( ( middle - Transpose(middle_rightfactor) * middle_rightfactor ).issmall() );
    
    return middle_rightfactor * multiplier; //TODO Probelesen
}

