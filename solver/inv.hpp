#ifndef INCLUDEGUARD_SOLVER_INVERSEOPERATOR
#define INCLUDEGUARD_SOLVER_INVERSEOPERATOR


#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include "../basic.hpp"

#include "crm.hpp"
#include "cgm.hpp"




class InverseOperator final
: public LinearOperator 
{

    public:

        explicit InverseOperator() = delete;
        explicit InverseOperator( const InverseOperator& ) = delete;
        explicit InverseOperator( InverseOperator&& ) = delete;
        InverseOperator& operator=( const InverseOperator& vec ) = delete;
        InverseOperator& operator=( InverseOperator&& vec ) = delete; 

        
        explicit InverseOperator( const LinearOperator& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), op( op ) { 
            assert( op.getdimin() == op.getdimout() );
            LOG << "Inverse created" << ""; 
        }
        
        virtual ~InverseOperator() { 
            LOG << "Inverse destroyed" << "";
        }

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<InverseOperator> cloned = std::make_shared<InverseOperator>( op );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<InverseOperator> heir = std::make_unique<InverseOperator>( op );
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
            
            {
                dest.zero();
                
                ConjugateResidualMethod Solver( op );
                
                Solver.max_iteration_count *= 4;
                Solver.print_modulo = Solver.max_iteration_count;
                
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
    
};
  
  
inline InverseOperator inv( const LinearOperator& op )
{
    op.check();
    
    return InverseOperator( op );
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
