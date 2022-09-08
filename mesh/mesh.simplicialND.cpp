
#include <algorithm>
#include <sstream>
#include <map>
#include <string>
#include <utility>
#include <vector>


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
    data_subsimplices( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    data_firstparents( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    data_nextparents ( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    
    flags_simplices( innerdim+1, std::vector<SimplexFlag>() )
{
    rebuild();
    MeshSimplicialND::check();
}


MeshSimplicialND::MeshSimplicialND( const Mesh& mesh )
:
    Mesh( mesh.getinnerdimension(), mesh.getouterdimension() ),
    
    counter_simplices( mesh.getinnerdimension()+1 ),
    data_subsimplices( (mesh.getinnerdimension()+2) * (mesh.getinnerdimension()+1) / 2, std::vector<int>() ),
    data_firstparents( (mesh.getinnerdimension()+2) * (mesh.getinnerdimension()+1) / 2, std::vector<int>() ),
    data_nextparents ( (mesh.getinnerdimension()+2) * (mesh.getinnerdimension()+1) / 2, std::vector<int>() ),
    
    flags_simplices( mesh.getinnerdimension()+1, std::vector<SimplexFlag>() )
{
    
    /* load coordinates */
    
    getcoordinates() = mesh.getcoordinates();
    
    
    /* check input data */
    
    assert( getinnerdimension()             == mesh.getinnerdimension() );
    assert( getouterdimension()             == mesh.getouterdimension() );
    assert( getcoordinates().getdimension() ==      getouterdimension() );
    assert( getcoordinates().getnumber()     >                        0 );
    
    
    /* use input data */
    
    const int innerdim = getinnerdimension();
    std::vector<int>& ref = data_subsimplices[index_from_pair(innerdim,0)];
    
    counter_simplices[innerdim] = mesh.count_simplices( innerdim );
    counter_simplices[0]        = mesh.count_simplices( 0 );
    
    ref.resize( mesh.count_simplices( innerdim ) * ( innerdim+1 ), nullindex );
    
    for( int S = 0; S <  mesh.count_simplices( innerdim ); S++ )
    for( int i = 0; i <=                         innerdim; i++ )
      ref[ S * ( innerdim + 1 ) + i ] 
        = mesh.get_subsimplex( innerdim, 0, S, i );
    
    
    
    /* check more stuff */
    
    assert( counter_simplices[0]-1 == *std::max_element( ref.begin(), ref.end() ) );
    
    assert( ref.size() > 0 );
    
    assert( ref.size() % (getinnerdimension()+1) == 0 );
    
    assert( getcoordinates().getnumber()-1 == *std::max_element( ref.begin(), ref.end() ) );
    
    
    /* rebuild and check */
    
    rebuild();
    
    MeshSimplicialND::check();
    
}


MeshSimplicialND::MeshSimplicialND( 
    int innerdim,
    int outerdim,
    const Coordinates& coords,
    const std::vector<int>& simplex_vertices
)
:
    Mesh( innerdim, outerdim ),
    
    counter_simplices( innerdim+1, 0 ),
    data_subsimplices( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    data_firstparents( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    data_nextparents ( (innerdim+2) * (innerdim+1) / 2, std::vector<int>() ),
    
    flags_simplices( innerdim+1, std::vector<SimplexFlag>() )   
{
    
    /* load coordinates */
    coords.check();
    getcoordinates() = coords;
    
    
    /* check input data */
    assert( innerdim >= 0 && outerdim >= 0 );
    assert( getcoordinates().getdimension() == outerdim );
    assert( getcoordinates().getdimension() == getouterdimension() );
    assert( getcoordinates().getnumber() > 0 );
    assert( simplex_vertices.size() > 0 );
    assert( getcoordinates().getnumber() - 1 == *std::max_element( simplex_vertices.begin(), simplex_vertices.end() ) );
    assert( simplex_vertices.size() % (innerdim+1) == 0 );
    
    
    /* use input data */
    counter_simplices[innerdim]                    = simplex_vertices.size() / (innerdim+1);
    counter_simplices[0]                           = getcoordinates().getnumber();
    data_subsimplices[index_from_pair(innerdim,0)] = simplex_vertices;
    
    
    /* check more stuff */
    assert( counter_simplices[0]-1 == *std::max_element( simplex_vertices.begin(), simplex_vertices.end() ) );
    
    
    LOG << "Rebuild inside constructor..." << space << innerdim << space << outerdim << space << nl;
    rebuild();
    LOG << "...done" << nl;
    
    /* check and exit */
    MeshSimplicialND::check();
    
}

MeshSimplicialND::~MeshSimplicialND()
{
    MeshSimplicialND::check();    
}



bool MeshSimplicialND::compare ( const MeshSimplicialND& mesh ) const 
{
  if( getinnerdimension() != mesh.getinnerdimension() ) return false;
  if( getouterdimension() != mesh.getouterdimension() ) return false;
  
  if( getcoordinates()    != mesh.getcoordinates()    ) return false;
  
  for( int d = 0; d <= getinnerdimension(); d++ )
    if( count_simplices(d) != mesh.count_simplices(d) ) return false;
  
  /* comparison of subsimplex lists and neighbor lists */
  
  for( int sup = 1; sup <= getinnerdimension(); sup++ )
  for( int sub = 0; sub <                  sup; sub++ )
    for( int S = 0; S < count_simplices(sup); S++ )
      for( int si = 0; si < count_subsimplices(sup,sub); si++ )
        if( data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ]
            !=
            mesh.data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ]
            ||
            data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ]
            !=
            mesh.data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ]
          )
          return false;
  
  /* comparison of first parent simplices */
  
  for( int sup = 1; sup <= getinnerdimension(); sup++ )
  for( int sub = 0; sub <                  sup; sub++ )
    for( int s = 0; s < count_simplices(sub); s++ )
      if( data_firstparents[index_from_pair(sup,sub)][ s ]
          !=
          mesh.data_firstparents[index_from_pair(sup,sub)][ s ]
        )
          return false;
  
  /* comparison of flags */
  
  for( int d = 0; d <= getinnerdimension(); d++ )
      for( int s = 0; s < count_simplices(d); s++ )
        if( flags_simplices[d][s] != mesh.flags_simplices[d][s] )
            return false;
          
          
  return true;
  
}
        
        

void MeshSimplicialND::check() const
{
    #ifdef NDEBUG
    return;
    #endif
    
    /* check that the number of arrays makes sense */
    
    assert( data_subsimplices.size() == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    assert( data_firstparents.size() == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    assert( data_nextparents.size()  == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    
    assert( flags_simplices.size() == getinnerdimension() + 1 );
    
    /* check that the dimensions of the arrays make sense */
    
    assert( counter_simplices.size() == getinnerdimension() + 1 );
    
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      assert( data_subsimplices[index_from_pair(sup,sub)].size() == counter_simplices[sup] * count_subsimplices(sup,sub) );
      assert( data_nextparents [index_from_pair(sup,sub)].size() == counter_simplices[sup] * count_subsimplices(sup,sub) );
      assert( data_firstparents[index_from_pair(sup,sub)].size() == counter_simplices[sub]                               );
    }
    
    for( int d = 0; d <= getinnerdimension(); d++ )
      assert( flags_simplices[d].size() == counter_simplices[d] );


    /* check the simplex counter function */
    
    for( int dim = 0; dim <= getinnerdimension(); dim++ )
      assert( counter_simplices[dim] == count_simplices(dim) );
    
    /* check the subsimplex lists:
     *   - non nullindex
     *   - array bounds
     *   - no duplicates
     */
    
    for( int sup = 1; sup <=         getinnerdimension(); sup++ )
    for( int sub = 0; sub <                          sup; sub++ )
    for( int S   = 0; S   <          count_simplices(sup); S++   )
    for( int si  = 0; si  <  count_subsimplices(sup,sub); si++  )
    {
      assert( data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] != nullindex          );
      assert( data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] >= 0                  );
      assert( data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] < count_simplices(sub) );
      
      for( int ti = 0; ti < si; ti++ )
        assert( data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] 
                !=
                data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + ti ]
              );
    }
    
    /* check the subsimplex lists:
     *   - transitivity of simplex 
     */
    
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,sub), IndexRange(0,sup) );
      
      for( int S  = 0; S  <         count_simplices(sup);  S++ )
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
    
    for( int sup = 1; sup <=         getinnerdimension(); sup++ )
    for( int sub = 0; sub <                          sup; sub++ )
    for( int S   = 0; S   <          count_simplices(sup); S++   )
    for( int si  = 0; si  <  count_subsimplices(sup,sub); si++  )
    {
      if( data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] == nullindex )
        continue;
      
      LOG << sup << space << sub << space << S << space << si << space 
                << data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] << nl;     
      assert( data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] >= 0                   );
      assert( data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ] <  count_simplices(sup) );
    }
    
    
    /* check the first parent lists:
     *   - check whether listed first parent is actually a parent
     *   - check successive parents
     */
    
    for( int sup = 1; sup <=   getinnerdimension(); sup++ )
    for( int sub = 0; sub <                    sup; sub++ )
    for( int s   = 0; s   < counter_simplices[sub]; s++   )
    {
      int S = data_firstparents[index_from_pair(sup,sub)][s];
      
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
    
    for( int sup = 1; sup <=         getinnerdimension(); sup++ )
    for( int sub = 0; sub <                          sup; sub++ )
    for( int S   = 0; S   <       counter_simplices[sup]; S++   )
    for( int si  = 0; si  <  count_subsimplices(sup,sub); si++  )
    {
      int s = data_subsimplices[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + si ];
      
      int P = data_firstparents[index_from_pair(sup,sub)][s];
      
      while( P != S )
      {
        assert( is_subsimplex( sup, sub, P, s ) );
        int pi = get_subsimplex_index( sup, sub, P, s );
        P = get_nextparent_by_localindex( sup, sub, P, pi );
      }
        
      assert( P == S );
    }
    
    Mesh::check();
    
}






std::string MeshSimplicialND::text() const
{
    std::ostringstream os;
    
    os << "Triangulation of N-dimensional Manifold!" << nl;
    
    os << "inner dimension: " << getinnerdimension() << nl
       << "outer dimension: " << getouterdimension() << nl;
    
    os << "counting simplices per dimension" << nl;
    for( int d = 0; d <= getinnerdimension(); d++ )
      os << d << tab << count_simplices(d) << nl;
      
    os << "subsimplices" << nl;
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      
      os << sup << " -> " << sub << nl;
      for( int S =  0; S < count_simplices(sup); S++ )
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
      for( int S =  0; S < count_simplices(sup); S++ )
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
      for( int s =  0; s < count_simplices(sub); s++ )
        os << get_firstparent_of_subsimplex(sup,sub,s) << space;
      os << nl;
      
    }
    
    getcoordinates().print( os );
    
    os << "Finished printing" << nl;

    return os.str();
}






bool MeshSimplicialND::dimension_counted( int dim ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    return true;
}

int MeshSimplicialND::count_simplices( int dim ) const
{
  assert( 0 <= dim && dim <= getinnerdimension() );
  return counter_simplices[dim];
}

bool MeshSimplicialND::subsimplices_listed( int sup, int sub ) const
{
  assert( 0 <= sub && sub <= sup && sup <= getinnerdimension() );
  return true;
}

IndexMap MeshSimplicialND::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub <= sup && sup <= getinnerdimension() );
  assert( 0 <= cell && cell < counter_simplices[sup] );
  
  const int C = count_subsimplices( sup, sub );
  assert( C > 0 );
  
  if( sup == sub )
    
    return IndexMap( IndexRange(0,0), IndexRange(0,counter_simplices[sup]-1), { cell } );
    
  else
    
    return IndexMap(
      IndexRange( 0, C-1                      ), 
      IndexRange( 0, counter_simplices[sub]-1 ), 
      std::vector<int>(
        data_subsimplices[ index_from_pair(sup,sub) ].begin() + C * ( cell + 0 ), 
        data_subsimplices[ index_from_pair(sup,sub) ].begin() + C * ( cell + 1 )
        )
      );
  
  unreachable();
    
}

bool MeshSimplicialND::supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= getinnerdimension() );
    return true;
}

const std::vector<int> MeshSimplicialND::getsupersimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub <= sup && sup <= getinnerdimension() );
  assert( 0 <= cell && cell < counter_simplices[sub] );
  
  if( sup == sub ) return { cell };
  
  std::vector<int> temp(0);
  int S = data_firstparents[index_from_pair(sup,sub)][cell];
  
  assert( S != nullindex );
  
  while( S != nullindex )
  {
    assert( is_subsimplex( sup, sub, S, cell ) );
    temp.push_back( S );
    int celli = get_subsimplex_index( sup, sub, S, cell );
    S = data_nextparents[index_from_pair(sup,sub)][ S * count_subsimplices(sup,sub) + celli ]; 
  }
  return temp;
}




SimplexFlag MeshSimplicialND::get_flag( int dim, int cell ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= cell && cell < count_simplices( dim ) );
    return flags_simplices[dim][cell];
}

void MeshSimplicialND::set_flag( int dim, int cell, SimplexFlag flag ) 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= cell && cell < count_simplices( dim ) );
    flags_simplices[dim][cell] = flag;
}





void MeshSimplicialND::rebuild()
{
    
    /* we need that the outer arrays have the right size */
    
    assert( data_subsimplices.size() == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    assert( data_firstparents.size() == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    assert( data_nextparents.size()  == (getinnerdimension()+2) * (getinnerdimension()+1) / 2 );
    
    
    /* first we clean all the arrays except for the volume->vertex array
       we also forget the old simplex counters except for the extremal dimensions */
    
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    {
      
      /* remove the simplex-subsimplex lists except for the original one */
    
      for( int sub = 0; sub <                  sup; sub++ )
      {
        if( sup == getinnerdimension() && sub == 0 ) continue;
        
        data_firstparents[ index_from_pair(sup,sub) ].resize(0);
        data_nextparents [ index_from_pair(sup,sub) ].resize(0);
        data_subsimplices[ index_from_pair(sup,sub) ].resize(0);
        
      }
    
      /* reset the simplex counter except for top-dimensional (and vertices) */
    
      if( sup < getinnerdimension() ) counter_simplices[sup] = 0;
    
    }
    
    
    /* rebuild the simplex->vertex arrays */
    
    for( int sub = 1; sub < getinnerdimension(); sub++ )
    {
      
      /* sigmas as auxiliary construction */
    
      std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,sub), IndexRange(0,getinnerdimension()) );
      
      /* collect the raw data */
      
      std::vector<std::vector<int>> temp( counter_simplices[getinnerdimension()] * sigmas.size(), 
                                          std::vector<int>(sub+1,nullindex)
                                        );
      
      for( int S  = 0; S  < counter_simplices[getinnerdimension()];  S++ )
      for( int si = 0; si <                          sigmas.size(); si++ )
      for( int vi = 0; vi <                                  sub+1; vi++ )
        temp.at( S * sigmas.size() + si ).at(vi) 
          = data_subsimplices[index_from_pair(getinnerdimension(),0)][ S * (getinnerdimension()+1) + sigmas[si][vi] ];
        
      std::sort( temp.begin(), temp.end() );
      auto it = std::unique( temp.begin(), temp.end() );
      temp.resize( it - temp.begin() );
      
      /* fill in the sorted data */
      
      data_subsimplices[index_from_pair(sub,0)].resize( temp.size() * (sub+1) );
      for( int s  = 0; s  <  temp.size();  s++ )
      for( int vi = 0; vi <=         sub; vi++ )
        data_subsimplices[index_from_pair(sub,0)].at( s * (sub+1) + vi ) = temp.at(s).at(vi);
      
      /* check those data */
      
      for( int i = 0; i < data_subsimplices[index_from_pair(sub,0)].size(); i++ ) {
        assert( data_subsimplices[index_from_pair(sub,0)].at(i) != nullindex            );
        assert( data_subsimplices[index_from_pair(sub,0)].at(i) >= 0                    );
        assert( data_subsimplices[index_from_pair(sub,0)].at(i) <  counter_simplices[0] );
      }
      
    }
    
    
    /*
     * reset the simplex counters 
     */
    for( int dim = 1; dim < getinnerdimension(); dim++ )
      counter_simplices[dim] = data_subsimplices[index_from_pair(dim,0)].size() / (dim+1);
     
    /*
     * resize arrays for next/parents for all. 
     * resize the subsimplex array for all except those to the vertices. 
     */
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 0; sub <                  sup; sub++ )
    {
      data_nextparents [index_from_pair(sup,sub)].resize( counter_simplices[sup] * count_subsimplices(sup,sub), nullindex );
      data_firstparents[index_from_pair(sup,sub)].resize( counter_simplices[sub],                               nullindex );
      
      if( sub == 0 ) continue;
      
      data_subsimplices[index_from_pair(sup,sub)].resize( counter_simplices[sup] * count_subsimplices(sup,sub), nullindex );  
    }
    
    
    /* fill up the firstparent and nextparent lists of vertices */
    
    for( int sup = 1; sup <=       getinnerdimension(); sup++ )
    for( int S   = 0; S   <     counter_simplices[sup];   S++ )
    for( int vi  = 0; vi  <=                       sup;  vi++ )
    {
      int v = data_subsimplices[index_from_pair(sup,0)][ S * (sup+1) + vi ];
      
      assert( v != nullindex );
      assert( 0 <= v && v < counter_simplices[0] );
      
      int old_first_parent = data_firstparents[index_from_pair(sup,0)][v];
      
      data_firstparents[index_from_pair(sup,0)][v] = S;
      
      assert( data_nextparents[index_from_pair(sup,0)][ S * (sup+1) + vi ] == nullindex );
      
      data_nextparents[index_from_pair(sup,0)][ S * (sup+1) + vi ] = old_first_parent;
    }
    
    
    /* fill up the firstparent and nextparent lists */
    /* fill also the subsimplex lists               */
    
    for( int sup = 1; sup <= getinnerdimension(); sup++ )
    for( int sub = 1; sub <                  sup; sub++ )
    {
      std::vector<IndexMap> sigmas = generateSigmas( IndexRange(0,sub), IndexRange(0,sup) );
      
      assert( sigmas.size() == count_subsimplices(sup,sub) );
      
      for( int S  = 0; S  < counter_simplices[sup]; S++ )
      for( int s  = 0; s  < counter_simplices[sub]; s++ )
      for( int si = 0; si < count_subsimplices(sup,sub); si++ )
      {
        std::vector<int> S_vertices( data_subsimplices[index_from_pair(sup,0)].begin() + (S  ) * (sup+1),
                                     data_subsimplices[index_from_pair(sup,0)].begin() + (S+1) * (sup+1)
                                   );
        
        std::vector<int> s_vertices( data_subsimplices[index_from_pair(sub,0)].begin() + (s  ) * (sub+1),
                                     data_subsimplices[index_from_pair(sub,0)].begin() + (s+1) * (sub+1)
                                   );
        
        IndexMap S_vertices_im( IndexRange(0,sup), IndexRange(0,counter_simplices[0]-1), S_vertices );
        IndexMap s_vertices_im( IndexRange(0,sub), IndexRange(0,counter_simplices[0]-1), s_vertices );
        
        if( s_vertices_im == S_vertices_im * sigmas[si] ){
          
          /* update the subsimplex list of the parent simplex */ 
          
          data_subsimplices[index_from_pair(sup,sub)].at( S * count_subsimplices(sup,sub) + si ) = s;
          
          /* update first parent of subsimplex and next parent of supersimplex */
          
          int old_first_parent = data_firstparents[ index_from_pair(sup,sub) ][ s ];
          data_firstparents[ index_from_pair(sup,sub) ][ s ] = S;
          assert( data_nextparents[ index_from_pair(sup,sub) ][ S * count_subsimplices(sup,sub) + si ] == nullindex );
          data_nextparents[ index_from_pair(sup,sub) ][ S * count_subsimplices(sup,sub) + si ] = old_first_parent;
          
        }
        
      }
      
    }
    
    // TODO: set flags!
    
}



void MeshSimplicialND::bisect_edge( int e )
{
    check(); assert( false );
}


void MeshSimplicialND::uniformrefinement()
{
    const int M = poweroftwo( getinnerdimension() );
    
    /*
     * iterate over the edges and create new vertices
     * integrate vertices into the coordinate block
     */
    
    getcoordinates().addcoordinates( counter_simplices[1] );
    
    for( int e = 0; e < counter_simplices[1]; e++ )
      getcoordinates().loadvector( 
        counter_simplices[0] + e, get_simplex_midpoint( 1, e )
      );
    
    /*
     * working copy of the array for top-dimensional simplices 
     * auxiliary save the pattern for new simplices 
     * iterate over main array to create the new top-dimensional simplices 
     */
    
    std::vector<int> temp( counter_simplices[getinnerdimension()] * M, nullindex );
    
    for( int S = 0; S < counter_simplices[getinnerdimension()]; S++ )
    for( int i = 0; i <                                      M; i++ )
    {
      
      std::vector<int> code1( getinnerdimension() + 1, 0 );
      std::vector<int> code2( getinnerdimension() + 1, 0 );
      for( int v = 0; v < getinnerdimension(); v++ )
      {
        // ( i >> v ) % 2; // TODO
      }
      
    }
    
    /*
     * Update counters for vertices and volumes 
     */
    
    counter_simplices[0] = counter_simplices[0] + counter_simplices[1];
    
    counter_simplices[ getinnerdimension() ] = M * counter_simplices[ getinnerdimension() ];
    
    /*
     * Rebuild 
     */
    
    rebuild();
    
    /*
     * Check 
     */
    
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






