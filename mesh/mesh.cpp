

#include <algorithm>


#include "../combinatorics/generateindexmaps.hpp"
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
  for( int sup = 0; sup <= innerdimension; sup++ )
  for( int sub = 0; sub <= sup; sub++ )
  {
//     std::cout << sup << space << sub << std::endl;
    IndexRange from( 0, countsubsimplices( sup, sub ) - 1 );
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





int Mesh::index_from_pair( int maxdimension, int sup, int sub ) 
{
  assert( 0 <= sub && sub < sup && sup <= maxdimension );
  return sub + (sup-1) * sup / 2;
}

void Mesh::index_to_pair( int maxdimension, int index, int& sup, int& sub ) 
{
  // TODO: improve this incredibly ineffecient computation.
  assert( 0 <= sub && sub <= sup && sup <= maxdimension );
  for( int _sup = 1; _sup <= maxdimension; _sup++ )
  for( int _sub = 0; _sub < sup; _sub++ )
    if( index == index_from_pair( maxdimension, _sup, _sub ) ){
      sup = _sup; sub = _sub; return;
    }
}

int Mesh::count_subsimplices( int maxdimension, int sup, int sub )
{
  assert( 0 <= sub && sub <= sup && sup <= maxdimension );
  return binomial( sup + 1, sub + 1);
}




bool Mesh::is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im.rangecontains( cellsub );
  
}

int  Mesh::get_subsimplix_index( int sup, int sub, int cellsup, int cellsub ) const
{
  
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im.rangeposition( cellsub );
  
}

int Mesh::get_subsimplix( int sup, int sub, int cellsup, int localindex ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im[ localindex ];  
}

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


int Mesh::get_supersimplex_index( int sup, int sub, int cellsup, int cellsub ) const
{
  
  assert( supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  return it - parents.begin();
  
}



