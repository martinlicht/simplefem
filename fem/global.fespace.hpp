#ifndef INCLUDEGUARD_FEM_GLOBAL_FESPACES
#define INCLUDEGUARD_FEM_GLOBAL_FESPACES




/* 
 * 
 * Global FE Space 
 * 
 * 
 * 
 * 
 */




class BrokenFESpace
{
  virtual SimplexFESpace get_local_FESpace( int t ) const = 0;
  
  // TODO: Submesh restriction 
  virtual BrokenFESpace get_derivative_space() const = 0;
  virtual BrokenFESpace get_jump_space() const = 0;
  
}

class UniformBrokenFESpace
: public virtual BrokenFESpace
{
  virtual SimplexFESpace get_local_FESpace( int t ) const;  
}

class GeneralBrokenFESpace
: public virtual BrokenFESpace
{
  virtual SimplexFESpace get_local_FESpace( int t ) const;
}


class ConformingFESpace
{
  
  virtual SimplexFESpace get_local_contribution( int dim, int simplex ) const = 0;
  virtual SimplexFESpace get_local_space( int dim, int simplex ) const = 0;
  
  // TODO: Submesh restriction
  virtual ConformingFESpace get_trace_space( int dim ) const = 0;
  virtual ConformingFESpace get_derivative_space( int dim ) const = 0;
  
}

class UniformConformingFESpace
: public virtual ConformingFESpace
{
  virtual SimplexFESpace get_local_contribution( int dim, int simplex ) const;
  virtual SimplexFESpace get_local_space( int dim, int simplex ) const;
}

class GeneralConformingFESpace
: public virtual ConformingFESpace
{
  virtual SimplexFESpace get_local_contribution( int dim, int simplex ) const;
  virtual SimplexFESpace get_local_space( int dim, int simplex ) const;
}





#endif