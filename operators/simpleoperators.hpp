#ifndef INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS_HPP
#define INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS_HPP

#include <cassert>
#include <functional>
#include <memory>
#include <ostream>
#include <utility>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Identity Operators 
****  
************************/


class IdentityOperator final
: public LinearOperator
{

    public:

        IdentityOperator()                                         = delete;
        IdentityOperator( const IdentityOperator& )                = default;
        IdentityOperator( IdentityOperator&& )                     = default;
        IdentityOperator& operator=( const IdentityOperator& vec ) = default;
        IdentityOperator& operator=( IdentityOperator&& vec )      = default; 
        
        explicit IdentityOperator( int n );
        virtual ~IdentityOperator();

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<IdentityOperator> cloned = std::make_shared<IdentityOperator>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<IdentityOperator> heir = std::make_unique<IdentityOperator>( *this );
            return heir;
        }
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual void print( std::ostream& ) const override;

        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
    
};
  
  
inline IdentityOperator operator*( const IdentityOperator& left, const IdentityOperator& right )
{
    left.check();
    right.check();
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    return IdentityOperator( left.getdimin() );
}  




/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class ScalingOperator final
: public LinearOperator /* every scaling operation is a linear operator */
{

    public:

        ScalingOperator()                                        = delete;
        ScalingOperator( const ScalingOperator& )                = default;
        ScalingOperator( ScalingOperator&& )                     = default;
        ScalingOperator& operator=( const ScalingOperator& vec ) = default;
        ScalingOperator& operator=( ScalingOperator&& vec )      = default; 
        
        explicit ScalingOperator( int, Float s );
        virtual ~ScalingOperator();

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<ScalingOperator> clone = std::make_shared<ScalingOperator>( *this );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<ScalingOperator> heir = std::make_unique<ScalingOperator>( *this );
            return heir;
        }
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual void print( std::ostream& ) const override;

        Float getscaling() const;
        void setscaling( Float s );

        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:

        Float scaling;
    
};
  
  
inline ScalingOperator operator*( const ScalingOperator& left, const ScalingOperator& right )
{
    left.check();
    right.check();
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    return ScalingOperator( left.getdimout(), left.getscaling() * right.getscaling() );
}  




/************************
****
****  Class for diagonal matrices 
****  - instantiates LinearOperator
****  - only square matrices 
****  
************************/

class DiagonalOperator final
: public LinearOperator 
{

    public:

        DiagonalOperator()                                         = delete;
        DiagonalOperator( const DiagonalOperator& )                = default;
        DiagonalOperator( DiagonalOperator&& )                     = default;
        DiagonalOperator& operator=( const DiagonalOperator& vec ) = default;
        DiagonalOperator& operator=( DiagonalOperator&& vec )      = default; 

        explicit DiagonalOperator( int, Float s );
        explicit DiagonalOperator( const FloatVector& dia );
        explicit DiagonalOperator( FloatVector&& dia );
        explicit DiagonalOperator( int, const ScalingOperator& scaling );
        explicit DiagonalOperator( int, const std::function<Float(int)>& );
        virtual ~DiagonalOperator();
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<DiagonalOperator> clone = std::make_shared<DiagonalOperator>( *this );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<DiagonalOperator> heir = std::make_unique<DiagonalOperator>( *this );
            return heir;
        }
        

        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual std::string text( const bool embellish ) const;
        
        virtual void print( std::ostream& ) const override;

        FloatVector& getdiagonal();
        const FloatVector& getdiagonal() const;
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
        
        const DiagonalOperator sqrt() const;

    private:

        FloatVector diagonal;
            
};
  
  

inline DiagonalOperator operator*( const DiagonalOperator& left, const DiagonalOperator& right )
{
    left.check();
    right.check();
    
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    const FloatVector& leftdia = left.getdiagonal();
    const FloatVector& rightdia = right.getdiagonal();
    const int dimension = leftdia.getdimension();
    
    assert( leftdia.getdimension() == rightdia.getdimension() );
    
    return DiagonalOperator( FloatVector( leftdia.getdimension(), 
                                          [&](int d) -> Float { 
                                            assert( 0 <= d && d < dimension ); 
                                            return leftdia[d] * rightdia[d];
                                          }
                                        )
                           );
}  
  
  
  
  

#endif
