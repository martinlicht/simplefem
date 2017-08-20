
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
    assert( data_simplex_subsimplices.size()               == getinnerdimension() * (getinnerdimension()+1) / 2 );
    assert( data_simplex_firstparent_simplex.size()        == getinnerdimension()                               );
    assert( data_simplex_nextparent_of_subsimplices.size() == getinnerdimension() * (getinnerdimension()+1) / 2 );
    
    /* check that the dimensions of the arrays make sense */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      assert( data_simplex_subsimplices[index_from_pair(sup,sub)].size() == counter_simplices[sup] );
      assert( data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)].size() == counter_simplices[sup] );
      assert( data_simplex_firstparent_simplex[index_from_pair(sup,sub)].size() == counter_simplices[sub] );
    }
    
    /* check the subsimplex lists:
     *   - non nullindex
     *   - array bounds
     *   - no duplicates
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
      for( int S = 0; S < countsimplices(sup); S++ )
        for( int si = 0; si < count_subsimplices(sup,sub); si++ )
        {
          assert( data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] != nullindex );
          assert( 0 <= data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] );
          assert( data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] < countsimplices(sub) );
          for( int ti = 0; ti < count_subsimplices(sup,sub); ti++ )
            if( si != ti )
              assert( data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] 
                      !=
                      data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + ti ]
                    );
        }
    
    /* check the subsimplex lists:
     *   - transitivity of simplex 
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,sub), IndexRange(0,sup) );
      
      for( int S = 0; S < countsimplices(sup); S++ )
      for( int si = 0; si < count_subsimplices(sup,sub); si++ )
      {
        int s = get_subsimplex( sup, sub, S, si );
        IndexMap sup_vertices = getsubsimplices( sup, 0, S );
        IndexMap sub_vertices = getsubsimplices( sub, 0, s );
        assert( sub_vertices == sup_vertices * sigmas[si] );
      }
    }
    
    
    
    
    /* check the neighborhood lists:
     *   - non nullindex
     *   - array bounds
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
      for( int S = 0; S < countsimplices(sup); S++ )
        for( int si = 0; si < count_subsimplices(sup,sub); si++ )
        {
          assert( data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] != nullindex );
          assert( 0 <= data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] );
          assert( data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] < countsimplices(sub) );
        }
    
    
    
    /* check the first parent lists:
     *   - check whether listed first parent is actually a parent
     *   - check successive parents
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
      for( int s = 0; s < counter_simplices[sub]; s++ )
      {
        int S = data_simplex_firstparent_simplex[index_from_pair(sup,sub)][s];
        
        assert( S != nullindex );
        
        while( S != nullindex )
        {
          assert( is_subsimplex( sup, sub, S, s ) );
          int si = get_subsimplex_index( sup, sub, S, s );
          S = get_nextparent_by_localindex( sup, sub, S, si );
        }
        
        assert( S == nullindex );
      }
    
    /* traverse the parent lists from each subsimplex:
     *   - check each one is actually listed as a parent
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    for( int S  = 0; S  <      counter_simplices[sup]; S++ )
    for( int si = 0; si < count_subsimplices(sup,sub); si++ )
    {
      int s = data_simplex_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ];
      
      int P = data_simplex_firstparent_simplex[index_from_pair(sup,sub)][s];
      
      while( P != S )
      {
        assert( is_subsimplex( sup, sub, P, s ) );
        int pi = get_subsimplex_index( sup, sub, P, s );
        S = get_nextparent_by_localindex( sup, sub, P, pi );
      }
        
      assert( P == S );
      
    }
    
    
    Mesh::check();
    
}






void MeshSimplicialND::print( std::ostream& os ) const
{
    os << "Printe Triangulation of N-dimensional Manifold!" << std::endl;
    
    os << "inner dimension: " << getinnerdimension() << nl
       << "outer dimension: " << getouterdimension() << nl;
    
    os << "counting simplices per dimension" << nl;
    for( int d = 0; d <= getinnerdimension(); d++ )
      os << d << tab << countsimplices(d) << nl;
      
    os << "subsimplices" << nl;
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      
      os << sup << " -> " << sub << nl;
      for( int S =  0; S < countsimplices(sup); S++ )
      {
        
        os << sup << ": ";
        for( int si = 0; si < count_subsimplices(sup,sub); si++ )
          os << get_subsimplex( sup, sub, S, si ) << space;
        os << nl;
        
      }
      
    }
    
    os << "next parents" << nl;
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      
      os << sup << " -> " << sub << nl;
      for( int S =  0; S < countsimplices(sup); S++ )
      {
        
        os << sup << ": ";
        for( int si = 0; si < count_subsimplices(sup,sub); si++ )
          os << get_nextparent_by_localindex( sup, sub, S, si ) << space;
        os << nl;
        
      }
      
    }
    
    os << "first parents" << nl;
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      
      os << sup << " -> " << sub << nl;
      for( int s =  0; s < countsimplices(sub); s++ )
        os << get_firstparent_of_subsimplex(sup,sub,s) << space;
      os << nl;
      
    }
    
    getcoordinates().print( os );
    
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

IndexMap MeshSimplicialND::getsubsimplices( int sup, int sub, int cell ) const
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






void MeshSimplicialND::rebuild()
{
    /* we need that the outer arrays have the right size */
    assert( data_simplex_subsimplices.size()               == getinnerdimension() * (getinnerdimension()+1) / 2 );
    assert( data_simplex_firstparent_simplex.size()        == getinnerdimension()                               );
    assert( data_simplex_nextparent_of_subsimplices.size() == getinnerdimension() * (getinnerdimension()+1) / 2 );
    
    /* first we clean all the arrays except for the volume->vertex array */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      if( sup == getinnerdimension() || sub == 0 ) continue;
      
      data_simplex_firstparent_simplex[index_from_pair(sup,sub)].resize(0);
      data_simplex_nextparent_of_subsimplices[index_from_pair(sup,sub)].resize(0);
      data_simplex_subsimplices[index_from_pair(sup,sub)].resize(0);
      
    }
    /* we also forget the old simplex counters except for the extremal dimensions */
    for( int dim = 1; dim < getinnerdimension(); dim++ )
      counter_simplices[dim] = 0;
    
    
    /* rebuild the simplex->vertex arrays */
    for( int sub = 0; sub < getinnerdimension(); sub++ )
    {
      std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,sub), IndexRange(0,getinnerdimension()) );
      
      std::vector<std::vector<int>> temp;
      temp.resize( counter_simplices[getinnerdimension() * sigmas.size() ], std::vector<int>(sub+1,nullindex) );
      
      for( int S  = 0; S  < counter_simplices[getinnerdimension()];  S++ )
      for( int si = 0; si <                          sigmas.size(); si++ )
      for( int vi = 0; vi <                                  sub+1; vi++ )
        temp[ S * sigmas.size() + si ][vi] = S * (getinnerdimension()+1) + sigmas[si][vi];
        
      std::sort( temp.begin(), temp.end() );
      auto it = std::unique( temp.begin(), temp.end() );
      temp.resize( it - temp.begin() );
      
      data_simplex_subsimplices[index_from_pair(sub,0)].resize( temp.size() * (sub+1) );
      for( int s  = 0; s  < temp.size();  s++ )
      for( int vi = 0; vi <=        sub; vi++ )
        data_simplex_subsimplices[index_from_pair(sub,0)][ s * (sub+1) + vi ] = temp[s][vi];
      
    }
    
    
    /* reset the simplex counters */
    for( int dim = 1; dim < getinnerdimension(); dim++ )
      counter_simplices[dim] = data_simplex_subsimplices[index_from_pair(dim,0)].size();
    
    
    /* reset the sizes of nextparent-lists */
    // TODO
    
    /* reset the sizes of firstparent-lists */
    
    
    /* fill up the firstparent and nextparent lists */
    
    
    check();
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






