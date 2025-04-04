#ifndef INCLUDEGUARD_SOLVER_INVERSEOPERATOR
#define INCLUDEGUARD_SOLVER_INVERSEOPERATOR


#include "../base/include.hpp"

#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "iterativesolver.hpp"
#include "sparsesolver.hpp"

#include "../sparse/matcsr.hpp"
#include "../sparse/rainbow.hpp"




template<typename Operator = LinearOperator>
class InverseOperator 
: public LinearOperator 
{

    public:

        InverseOperator()                                          = delete;
        InverseOperator( const InverseOperator& )                  = default;
        InverseOperator& operator=( const InverseOperator& invop ) = default;
        InverseOperator( InverseOperator&& )                       = default;
        InverseOperator& operator=( InverseOperator&& invop )      noexcept = default; 

        
        explicit InverseOperator( const Operator& op, Float precision, int print_modulo = -1 )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
          op( op ), precision(precision), print_modulo( print_modulo ), previous_sol( op.getdimin(), 0. )
        { 
            assert( op.getdimin() == op.getdimout() );
            
            if( print_modulo >= 0 ) {
                LOG << "Inverse created" << nl; 
            } 
        }
        
        virtual ~InverseOperator() noexcept = default;

        virtual InverseOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }

        virtual void check() const override { 
            op.check();
            previous_sol.check();
            assert( op.getdimin() == op.getdimout() );
            assert( op.getdimin() == previous_sol.getdimension() );
        }
        
        virtual std::string text() const override { 
            return "Inverse Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n" 
                    + tab_each_line( op.text() );
        }
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const override;
        
    protected:

        const Operator& op;
        Float precision;
        int print_modulo;

        mutable FloatVector previous_sol;

    public:
        bool use_previous_sol = true;
    
};

class PseudoInverseOperator
: public InverseOperator<LinearOperator> 
{

    public:

        PseudoInverseOperator()                                                = delete;
        PseudoInverseOperator( const PseudoInverseOperator& )                  = default;
        PseudoInverseOperator& operator=( const PseudoInverseOperator& invop ) = default;
        PseudoInverseOperator( PseudoInverseOperator&& )                       = default;
        PseudoInverseOperator& operator=( PseudoInverseOperator&& invop )      noexcept = default; 

        
        explicit PseudoInverseOperator( const LinearOperator& op, Float precision, int print_modulo = -1 )
        : InverseOperator( op, precision, print_modulo )
        { 
            assert( op.getdimin() == op.getdimout() );

            if( print_modulo >= 0 ) {
                LOG << "PseudoInverse created" << nl; 
            }
        }
        
        virtual ~PseudoInverseOperator() noexcept = default;

        virtual PseudoInverseOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }

        virtual std::string text() const override { 
            return "PseudoInverse Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n" 
                    + tab_each_line( op.text() );
        }
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const override;
        
    
};



template<typename T>
void InverseOperator<T>::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    if( use_previous_sol ) 
        dest = previous_sol;
    else 
        dest.zero();

    ConjugateResidualMethod solver( op );
    
    solver.max_iteration_count = op.getdimin();
    solver.print_modulo        = print_modulo;
    solver.verbosity           = ConjugateResidualMethod::VerbosityLevel::silent;
    
    solver.solve( dest, src );
    
    dest *= scaling;
    
    if( use_previous_sol ) previous_sol = dest;
    
}

template<>
void InverseOperator<MatrixCSR>::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    if( use_previous_sol ) 
        dest = previous_sol;
    else 
        dest.zero();

    const auto diagonal = op.getDiagonal();

    FloatVector res( dest );
    
    #if defined(_OPENMP)
    Rainbow rainbow( op );

    ConjugateGradientSolverCSR_Eisenstat_Rainbow( 
        src.getdimension(),
        dest.raw(), 
        src.raw(), 
        op.getA(), op.getC(), op.getV(), 
        res.raw(),
        precision,
        print_modulo,
        diagonal.raw(),
        1.0,
        rainbow.num_colors, rainbow.F.data(), rainbow.B.data(), rainbow.R.data()
    );
    #else
    ConjugateGradientSolverCSR_SSOR_Eisenstat( 
        src.getdimension(),
        dest.raw(), 
        src.raw(), 
        op.getA(), op.getC(), op.getV(), 
        res.raw(),
        precision,
        print_modulo,
        diagonal.raw(),
        1.0
    );
    #endif
    dest *= scaling;
    
    if( use_previous_sol ) previous_sol = dest;
    
}


void PseudoInverseOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    if( use_previous_sol ) 
        dest = previous_sol;
    else 
        dest.zero();

    ConjugateResidualMethod solver( op );
    
    solver.max_iteration_count = op.getdimin();
    solver.print_modulo        = print_modulo;
    solver.verbosity           = ConjugateResidualMethod::VerbosityLevel::silent;
    
    solver.solve( dest, src );
    
    dest *= scaling;
    
    if( use_previous_sol ) previous_sol = dest;
    
}


inline InverseOperator<LinearOperator> inv( const LinearOperator& op, Float precision, int print_modulo = -1 )
{
    op.check();
    return InverseOperator<LinearOperator>( op, precision, print_modulo );
} 
  
  
inline InverseOperator<MatrixCSR> inv( const MatrixCSR& op, Float precision, int print_modulo = -1 )
{
    op.check();
    return InverseOperator<MatrixCSR>( op, precision, print_modulo );
}  


inline PseudoInverseOperator pinv( const LinearOperator& op, Float precision, int print_modulo = -1 )
{
    op.check();
    return PseudoInverseOperator( op, precision, print_modulo );
} 
  
  




#endif
