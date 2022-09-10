#ifndef INCLUDEGUARD_SOLVER_ITERATIVESOLVER
#define INCLUDEGUARD_SOLVER_ITERATIVESOLVER


#include <ostream>
#include <limits>

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"


/************************
****
****  Abstract class for iterative solvers  
****  - uses iteration counter, error threshold, and internal residual vector 
****  
************************/


  
struct IterativeSolver
{
    
    enum class VerbosityLevel {
        silent = 0,
        resultonly = 1,
        verbose = 2
    };
    
    explicit IterativeSolver( const LinearOperator& A, Float threshold = desired_precision, int max_iteration_count = 0, int print_modulo = -1 )
    : A(A), 
      threshold( threshold ), 
      recent_deviation( 0. ), 
      max_iteration_count( max_iteration_count ),
      recent_iteration_count(0),
      print_modulo( print_modulo ), 
      verbosity( VerbosityLevel::verbose ) 
    {
        IterativeSolver::check();
    }
    
    virtual ~IterativeSolver() = default;

    virtual void check() const
    {
        assert( std::isfinite( threshold ) && threshold >= 0. );
        assert( std::isfinite( recent_deviation ) && recent_deviation >= 0. );
        
        assert( max_iteration_count >= 0 );
        assert( recent_iteration_count >= 0 );
        assert( recent_iteration_count <= max_iteration_count );

        assert( print_modulo >= -1 );
        
        A.check();
        assert( A.getdimin() == A.getdimout() );
    }

    virtual std::string text() const = 0;

    virtual void solve( FloatVector& unknown, const FloatVector& rhs ) const = 0;


    const LinearOperator& A;   

    mutable Float threshold;
    mutable Float recent_deviation;
    
    mutable int max_iteration_count;
    mutable int recent_iteration_count;
    
    mutable int print_modulo;
    
    mutable VerbosityLevel verbosity;
    
};



inline std::ostream& operator<<( std::ostream& os, const IterativeSolver& solver )
{
    os << solver.text();
    return os;
}














struct ConjugateGradientMethod
: public IterativeSolver
{

        using IterativeSolver::IterativeSolver;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;

};


struct ConjugateResidualMethod
: public IterativeSolver
{

        using IterativeSolver::IterativeSolver;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;

        virtual void solve_explicit( FloatVector&, const FloatVector& ) const;
        virtual void solve_robust( FloatVector&, const FloatVector& ) const;
        virtual void solve_fast( FloatVector&, const FloatVector& ) const;

};




struct PreconditionedConjugateResidualMethod
: public IterativeSolver
{

        explicit PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M );
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;
        
    private:

        const LinearOperator& M;

};


struct MinimumResidualMethod
: public IterativeSolver
{

        using IterativeSolver::IterativeSolver;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;

};

struct ResidualDescentMethod
: public IterativeSolver
{

        using IterativeSolver::IterativeSolver;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;

};

struct HerzogSoodhalterMethod
: public IterativeSolver
{

        using IterativeSolver::IterativeSolver;
        
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;

};


















 
  
  
  
  
  
#endif
