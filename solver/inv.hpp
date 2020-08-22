#ifndef INCLUDEGUARD_SOLVER_INVERSEOPERATOR
#define INCLUDEGUARD_SOLVER_INVERSEOPERATOR


#include <cassert>
#include <memory>
#include <ostream>


#include "../basic.hpp"

#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "crm.hpp"
#include "cgm.hpp"
#include "sparsesolver.hpp"

#include "../sparse/matcsr.hpp"




class InverseOperator final
: public LinearOperator 
{

    public:

        explicit InverseOperator() = delete;
        explicit InverseOperator( const InverseOperator& ) = delete;
        explicit InverseOperator( InverseOperator&& ) = delete;
        InverseOperator& operator=( const InverseOperator& vec ) = delete;
        InverseOperator& operator=( InverseOperator&& vec ) = delete; 

        
        explicit InverseOperator( const LinearOperator& op, Float tolerance )
        : LinearOperator( op.getdimout(), op.getdimin() ), op( op ), tolerance(tolerance) { 
            assert( op.getdimin() == op.getdimout() );
            LOG << "Inverse created" << ""; 
        }
        
        virtual ~InverseOperator() { 
            LOG << "Inverse destroyed" << "";
        }

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<InverseOperator> cloned = std::make_shared<InverseOperator>( op, tolerance );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<InverseOperator> heir = std::make_unique<InverseOperator>( op, tolerance );
            return heir;
        }

        virtual void check() const override { 
            op.check();
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Inverse Operator" << std::endl;
        }
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
            check();
            src.check();
            dest.check();
            
            assert( getdimin() == src.getdimension() );
            assert( getdimout() == dest.getdimension() );
            
            if( dynamic_cast<const MatrixCSR*>(&op) != nullptr ){
                
                dest.zero();
                
                FloatVector res( dest );
                
                auto opcsr = dynamic_cast<const MatrixCSR*>(&op);
                
                ConjugateGradientSolverCSR( 
                    src.getdimension(),
                    dest.raw(), 
                    src.raw(), 
                    opcsr->getA(), opcsr->getC(), opcsr->getV(), 
                    res.raw(),
                    tolerance,
                    src.getdimension()
                );
                
                dest *= scaling;
            
            } else {
            
                dest.zero();
                
                ConjugateResidualMethod Solver( op );
                
                Solver.max_iteration_count = op.getdimin();
                Solver.print_modulo        = Solver.max_iteration_count;
                Solver.verbosity           = ConjugateResidualMethod::VerbosityLevel::resultonly;
                
                Solver.solve_robust( dest, src );
                
                dest *= scaling;
            }
            
//             LOG << "call inverse";
//             op.print( std::cout );
//             op.apply( dest, src, scaling );    
//             LOG << "done inverse";
        }

    private:

        const LinearOperator& op;
        Float tolerance;
    
};
  
  
inline InverseOperator inv( const LinearOperator& op, Float tolerance )
{
    op.check();
    
    return InverseOperator( op, tolerance );
}  






//         using LinearOperator::apply;
//         void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
//             
//             dest.zero();
//             
//             ConjugateResidualMethod Solver( *op );
//             
//             Solver.max_iteration_count *= 4;
//             Solver.print_modulo = Solver.max_iteration_count;
//             
//             Solver.solve_robust( dest, src );
//             
//             dest *= scaling;
//         }








#endif
