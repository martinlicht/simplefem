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
    
    explicit IterativeSolver( Float threshold = desired_precision, int max_iteration_count = 0, int print_modulo = -1 )
    : threshold( threshold ), 
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

//         assert( print_modulo >= 0 );
    }

    virtual std::string text() const = 0;

//     // void lg() const { LOG << *this << std::endl; };

    virtual void solve( FloatVector& unknown, const FloatVector& rhs ) const = 0;

    mutable Float threshold;
    mutable Float recent_deviation;
    
    mutable int   max_iteration_count;
    mutable int   recent_iteration_count;
    
    mutable int   print_modulo;
    
    mutable VerbosityLevel verbosity;
    

};



inline std::ostream& operator<<( std::ostream& os, const IterativeSolver& solver )
{
    os << solver.text();
    return os;
}













/************************
****
****  Class for Conjugate Gradient Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class ConjugateGradientMethod
: public IterativeSolver
{

        public:
        
                explicit ConjugateGradientMethod( const LinearOperator& op );
                virtual ~ConjugateGradientMethod();

                virtual void check() const override;
                virtual std::string text() const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};
















/************************
****
****  Class for Conjugate Residual Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class ConjugateResidualMethod
: public IterativeSolver
{

        public:
        
                explicit ConjugateResidualMethod( const LinearOperator& op );
                virtual ~ConjugateResidualMethod();

                virtual void check() const override;
                virtual std::string text() const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

                virtual void solve_explicit( FloatVector&, const FloatVector& ) const;
                virtual void solve_robust( FloatVector&, const FloatVector& ) const;
                virtual void solve_fast( FloatVector&, const FloatVector& ) const;

        private: 

                const LinearOperator& A;   
};






















/************************
****
****  Class for Conjugate Residual Method
****  - instantiates IterativeSolver
****  - features iteration start and iteration step,
****    which can be called as such, or from solve().
****  
************************/




class PreconditionedConjugateResidualMethod
: public IterativeSolver
{

    public:
        
        explicit PreconditionedConjugateResidualMethod( const LinearOperator& op, const LinearOperator& M );
        virtual ~PreconditionedConjugateResidualMethod();

        virtual void check() const override;
        virtual std::string text() const override;
        
        virtual void solve( FloatVector&, const FloatVector& ) const override;
        
    private:

        const LinearOperator& A;   
        const LinearOperator& M;

};
  






















class MinimumResidualMethod
: public IterativeSolver
{

        public:
        
                explicit MinimumResidualMethod( const LinearOperator& op );
                virtual ~MinimumResidualMethod();

                virtual void check() const override;
                virtual std::string text() const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};















class ResidualDescentMethod
: public IterativeSolver
{

        public:
        
                explicit ResidualDescentMethod( const LinearOperator& op );
                virtual ~ResidualDescentMethod();

                virtual void check() const override;
                virtual std::string text() const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};























class HerzogSoodhalterMethod
: public IterativeSolver
{

        public:
        
                explicit HerzogSoodhalterMethod( const LinearOperator& op );
                virtual ~HerzogSoodhalterMethod();

                virtual void check() const override;
                virtual std::string text() const override;
                
                virtual void solve( FloatVector&, const FloatVector& ) const override;

        private: 

                const LinearOperator& A;   
};


















 
  
  
  
  
  
#endif
