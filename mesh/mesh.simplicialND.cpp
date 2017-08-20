
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"

#include "mesh.simplicialND.hpp"




MeshSimplicialND::MeshSimplicialND( int innerdim, int outerdim )
:
    Mesh( innerdim, outerdim ),
    
    counter_simplices( innerdim+1, 0 ),
    data_simplex_subsimplices( innerdim * (innerdim+1) / 2, std::vector<int>() ),
    data_simplex_firstparent_simplex( innerdim, std::vector<int>() ),
    data_simplex_nextparent_of_subsimplices( innerdim * (innerdim+1) / 2, std::vector<int>() )   
{
    rebuild();
    check();
}


MeshSimplicialND::MeshSimplicialND( 
    int innerdim,
    int outerdim,
    const Coordinates& coords,
    const std::vector<int> simplex_vertices
)
:
    Mesh( innerdim, outerdim ),
    
    counter_simplices( innerdim+1, 0 ),
    data_simplex_subsimplices( innerdim * (innerdim+1) / 2, std::vector<int>() ),
    data_simplex_firstparent_simplex( innerdim, std::vector<int>() ),
    data_simplex_nextparent_of_subsimplices( innerdim * (innerdim+1) / 2, std::vector<int>() )   
{
    
    getcoordinates() = coords;
    
    counter_simplices[innerdim] = simplex_vertices.size() / innerdim;
    
    rebuild();
    
    check();
}

MeshSimplicialND::~MeshSimplicialND()
{
    
}



bool MeshSimplicialND::operator== ( const MeshSimplicialND& mesh ) const 
{
  if( getinnerdimension() != mesh.getinnerdimension() ) return false;
  if( getouterdimension() != mesh.getouterdimension() ) return false;
  if( getcoordinates()    != mesh.getcoordinates()    ) return false;
  
  for( int d = 0; d <= getinnerdimension(); d++ )
    if( countsimplices(d) != mesh.countsimplices(d) ) return false;
  
  /* comparison of subsimplex lists and neighbor lists */
  
  for( int sup = 1; sup <= getinnerdimension(); sup++ )
  for( int sub = 0; sub <                  sup; sub++ )
    for( int S = 0; S < countsimplices(sup); S++ )
      for( int s = 0; s < count_subsimplices(sup,sub); s++ )
        if( data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + s ]
            !=
            mesh.data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + s ]
            ||
            data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + s ]
            !=
            mesh.data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + s ]
          )
          return false;
  
  /* comparison of first parent simplices */
  
  for( int sup = 1; sup <= getinnerdimension(); sup++ )
  for( int sub = 0; sub <                  sup; sub++ )
    for( int s = 0; s < countsimplices(sub); s++ )
      if( data_simplex_firstparent_simplex[index_from_pair(sup,sub)][ s ]
          !=
          mesh.data_simplex_firstparent_simplex[index_from_pair(sup,sub)][ s ]
        )
          return false;
  
  return true;
  
}
        
        

void MeshSimplicialND::check() const
{
    /* check that the number of arrays makes sense */
    
    /* check that the dimensions of the arrays make sense */
    
    /* check the subsimplex lists:
     *   - non nullindex
     *   - array bounds
     *   - no duplicates
     */
    
    /* check the neighborhood lists:
     *   - non nullindex
     *   - array bounds
     *   - no duplicates
     */
    
    /* check the first parent lists:
     *   - non nullindex
     *   - array bounds
     */
    
    /* traverse the parent lists from each subsimplex:
     *   - check each one is actually a parent
     */
    
    /* check the parent listentries:
     *   - check each one is actually a parent
     *   - check next parent 
     */
    
    Mesh::check();
    
}






void MeshSimplicialND::print( std::ostream& os ) const
{
    os << "Printe Triangulation of ND Manifold!" << std::endl;
    
    /* TODO: Printing */
    
    os << "Finished printing" << nl;
    
}






bool MeshSimplicialND::dimensioncounted( int dim ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    return true;
}

int MeshSimplicialND::countsimplices( int dim ) const
{
  assert( 0 <= dim && dim <= getinnerdimension() );
  return counter_simplices[dim];
}

bool MeshSimplicialND::subsimplices_listed( int sup, int sub ) const
{
  assert( 0 <= sub && sub < sup && sup <= getinnerdimension() );
  return true;
}

const IndexMap MeshSimplicialND::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub < sup && sup <= getinnerdimension() );
  assert( 0 <= cell && cell < counter_simplices[sup] );
  const int C = count_subsimplices( sup, sub );
  return IndexMap( IndexRange(0,1), IndexRange( 0, counter_simplices[sub]-1 ), 
                   std::vector<int>( data_simplex_subsimplices[ index_from_pair(sup,sub) ].begin() + cell * C, 
                                     data_simplex_subsimplices[ index_from_pair(sup,sub) ].begin() + cell * (C+1)
                                   )
                 );
    
}

bool MeshSimplicialND::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub < sup && sup <= getinnerdimension() );
    return true;
}

const std::vector<int> MeshSimplicialND::getsupersimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub < sup && sup <= getinnerdimension() );
  assert( 0 <= cell && cell < counter_simplices[sub] );
  // TODO: Write this method 
  return std::vector<int>( );
}










void MeshSimplicialND::bisect_edge( int e )
{
    check(); assert( e != e );
}


void MeshSimplicialND::uniformrefinement()
{
    check();
}








FloatVector MeshSimplicialND::get_simplex_midpoint( int dim, int s ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    FloatVector mid( getouterdimension() );
    
    for( int d = 0; d < getouterdimension(); d++ ) 
      mid[d] = 0.;
      
    for( int d = 0; d <  getouterdimension(); d++ )
    for( int c = 0; c <= getinnerdimension(); c++ )
      mid[d] += getcoordinates().getdata( get_subsimplex( dim, 0, s, c ), d );
    
    for( int d = 0; d < getouterdimension(); d++ ) 
      mid[d] /= getinnerdimension();
    
    return mid;
}






