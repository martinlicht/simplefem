

#include <algorithm>


#include "../combinatorics/generateindexmaps.hpp"
#include "mesh.hpp"
#include "../basic.hpp"


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
  for( int sup = 0; sup <= innerdimension; sup++ )
  for( int sub = 0; sub <= sup; sub++ )
  {
//     std::cout << sup << space << sub << std::endl;
    IndexRange from( 0, count_subsimplices( sup, sub ) - 1 );
    IndexRange to( 0, innerdimension );
//     std::cout << "generate sigmas" << nl << from << nl << to << std::endl;
    std::vector<IndexMap> sigmas = generateSigmas( from, to );
//     std::cout << "insert sigmas" << std::endl;
    auxdata[ std::pair<int,int>(sup,sub) ] = sigmas;
  }
  
//   std::cout << "mesh constructor" << std::endl;
  
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
  
  // check dimension 
  assert( outerdimension == coordinates.getdimension() );
  
  // the vertices and volumes must be counted 
  assert( dimensioncounted( innerdimension ) );
  assert( dimensioncounted( 0 ) );
  
  // the counting of the vertices must agree 
  assert( countsimplices(0) == coordinates.getnumber() );
  
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
  for( int _sub = 0; _sub <= sup; _sub++ )
    if( index == index_from_pair( _sup, _sub ) ){
      sup = _sup; sub = _sub; return;
    }
}

int Mesh::count_subsimplices( int sup, int sub ) const 
{
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  return binomial<int>( sup + 1, sub + 1);
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
  return im.rangeposition( cellsub );
  
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

