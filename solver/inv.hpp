#ifndef INCLUDEGUARD_SOLVER_INVERSEOPERATOR
#define INCLUDEGUARD_SOLVER_INVERSEOPERATOR


#include <cassert>
//#include <memory>


#include "../basic.hpp"

#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "iterativesolver.hpp"
#include "sparsesolver.hpp"

#include "../sparse/matcsr.hpp"




class InverseOperator final
: public LinearOperator 
{

    public:

        InverseOperator()                                          = delete;
        InverseOperator( const InverseOperator& )                  = default;
        InverseOperator& operator=( const InverseOperator& invop ) = default;
        InverseOperator( InverseOperator&& )                       = default;
        InverseOperator& operator=( InverseOperator&& invop )      = default; 

        
        explicit InverseOperator( const LinearOperator& op, Float tolerance, int print_modulo = -1 )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
          op( op ), tolerance(tolerance), print_modulo( print_modulo ), previous_sol( op.getdimin(), 0. )
        { 
            assert( op.getdimin() == op.getdimout() );
            
            LOG << "Inverse created" << nl; 
        }
        
        virtual ~InverseOperator() = default;

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
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
            check();
            src.check();
            dest.check();
            
            assert( getdimin() == src.getdimension() );
            assert( getdimout() == dest.getdimension() );
            
            if( use_previous_sol ) 
                dest = previous_sol;
            else 
                dest.zero();

            if( dynamic_cast<const MatrixCSR*>(&op) != nullptr ){
                
                const auto* opcsr = dynamic_cast<const MatrixCSR*>(&op);
                
                const auto diagonal = opcsr->diagonal();

                FloatVector res( dest );
                
                ConjugateGradientSolverCSR_SSOR( 
                    src.getdimension(),
                    dest.raw(), 
                    src.raw(), 
                    opcsr->getA(), opcsr->getC(), opcsr->getV(), 
                    res.raw(),
                    tolerance,
                    print_modulo,
                    diagonal.raw(),
                    1.0
                );
                
            } else {
            
                ConjugateGradientMethod Solver( op );
                
                Solver.max_iteration_count = op.getdimin();
                Solver.print_modulo        = print_modulo;
                Solver.verbosity           = ConjugateResidualMethod::VerbosityLevel::silent;
                
                Solver.solve( dest, src );
                
            }
            
            dest *= scaling;
            
            if( use_previous_sol ) 
                previous_sol = dest;
            
//             LOG << "call inverse" << nl;
//             LOG << op.text() << nl;   //op.print( std::cout );
//             op.apply( dest, src, scaling );    
//             LOG << "done inverse" << nl;
        }

    private:

        const LinearOperator& op;
        Float tolerance;
        int print_modulo;

        mutable FloatVector previous_sol;

    public:
        bool use_previous_sol = true;
    
};
  
  
inline InverseOperator inv( const LinearOperator& op, Float tolerance, int print_modulo = -1 )
{
    op.check();
    return InverseOperator( op, tolerance, print_modulo );
}  





#endif
