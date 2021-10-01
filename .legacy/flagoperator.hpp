#ifndef INCLUDEGUARD_OPERATOR_FLAGOPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_FLAGOPERATOR_HPP

#include <memory>
#include <vector>

#include "../basic.hpp"
#include "linearoperator.hpp"






/************************
****
****  Class for flagged operators 
****  - instantiates LinearOperator
****  
************************/

class FlagOperator final
: public LinearOperator 
{

    public:

        /* Constructors */
        
        explicit FlagOperator( const LinearOperator& op, const std::vector<bool>& destflag, const std::vector<bool>& srcflag );
        explicit FlagOperator( const LinearOperator& op, const std::vector<bool>& flag );
        
        /* standard interface */
        
        explicit FlagOperator() = delete;
        explicit FlagOperator( const FlagOperator& ) = default;
        explicit FlagOperator( FlagOperator&& ) = default;
        FlagOperator& operator=( const FlagOperator& vec ) = delete;
        FlagOperator& operator=( FlagOperator&& vec ) = delete;
        virtual ~FlagOperator();
        // TODO: Instantiate move semantics 

        /* standard methods for operators */
        
        virtual void check() const override;
        
        std::string text() const override;

        /* OTHER METHODS */
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<FlagOperator> cloned = std::make_shared<FlagOperator>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<FlagOperator> heir = std::make_unique<FlagOperator>( *this );
            return heir;
        }
        

        
        std::vector<bool>& getsrcflag();
        const std::vector<bool>& getsrcflag() const;
        
        std::vector<bool>& getdestflag();
        const std::vector<bool>& getdestflag() const;
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:

        const LinearOperator& op;

        std::vector<bool> destflag;
        std::vector<bool> srcflag;
            
};
  
  

  
  
  
  

#endif
