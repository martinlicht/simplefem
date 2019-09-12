#ifndef INCLUDEGUARD_FEM_ELEMENT_FESPACES
#define INCLUDEGUARD_FEM_ELEMENT_FESPACES


/* 
 * 
 * Simplex FE Space 
 * 
 */

struct SimplexFESpace
{
  
  int  dimension;
  int  polydegree;
  int  formdegree;
  bool trimmed;
  bool mathring;
  
  simplicialfemspace( int dim, int polydegree, int formdegree, bool trimmed, bool mathring )
  : dimension(dim), polydegree(polydegree), trimmed(trimmed), mathring(mathring)
  { };
  
  bool is_nullspace() const { unreachable() };
  
  bool is_scalar() const { unreachable() };
  
  int getdimension() const { return 0; };
  
}


inline is_comparable( SimplexFESpace fes1, SimplexFESpace fes2 )
{
  return ( fes1.dimension == fes2.dimension ) && ( fes1.formdegree == fes2.formdegree );
}


inline operator<( SimplexFESpace fes1, SimplexFESpace fes2 )
{
  assert( is_comparable( fes1, fes2 ) );
  // TODO: vergleich implementieren
}

inline operator>( SimplexFESpace fes1, SimplexFESpace fes2 )
{
  assert( is_comparable( fes1, fes2 ) );
  return fes2 < fes1;
}

inline operator<=( SimplexFESpace fes1, SimplexFESpace fes2 )
{
  assert( is_comparable( fes1, fes2 ) );
  return fes1 < fes2 && fes1 == fes2;
}

inline operator>=( SimplexFESpace fes1, SimplexFESpace fes2 )
{
  assert( is_comparable( fes1, fes2 ) );
  return fes2 < fes1 && fes1 == fes2;
}

/* TODO: Comparison operators for SimplexFESpaces */

/* TODO: Maximum/Minimum operators for SimplexFESpaces */







#endif