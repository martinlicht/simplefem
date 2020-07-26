#ifndef INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS
#define INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Identity Operators 
****  
************************/


class IdentityOperator:
public LinearOperator
{

    public:

        explicit IdentityOperator() = delete;
        explicit IdentityOperator( const IdentityOperator& ) = default;
        explicit IdentityOperator( IdentityOperator&& ) = default;
        IdentityOperator& operator=( const IdentityOperator& vec ) = default;
        IdentityOperator& operator=( IdentityOperator&& vec ) = default; 
        
        explicit IdentityOperator( int n ) : LinearOperator(n,n) { IdentityOperator::check(); };
        virtual ~IdentityOperator() { IdentityOperator::check(); };

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<IdentityOperator> clone = std::make_shared<IdentityOperator>( *this );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<IdentityOperator> heir = std::make_unique<IdentityOperator>( *this );
            return heir;
        }
        
        virtual void check() const override { LinearOperator::check(); };
        virtual void print( std::ostream& ) const override { LOG << "Print Identity Operator" << std::endl; };

        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override { dest = scaling * src; };
    
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


class ScalingOperator:
public LinearOperator /* every scaling operation is a linear operator */
{

    public:

        explicit ScalingOperator() = delete;
        explicit ScalingOperator( const ScalingOperator& ) = default;
        explicit ScalingOperator( ScalingOperator&& ) = default;
        ScalingOperator& operator=( const ScalingOperator& vec ) = default;
        ScalingOperator& operator=( ScalingOperator&& vec ) = default; 
        
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
        virtual void print( std::ostream& ) const override;

        Float getscaling() const;
        void setscaling( Float s );

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

class DiagonalOperator:
public LinearOperator 
{

    public:

        explicit DiagonalOperator() = delete;
        explicit DiagonalOperator( const DiagonalOperator& ) = default;
        explicit DiagonalOperator( DiagonalOperator&& ) = default;
        DiagonalOperator& operator=( const DiagonalOperator& vec ) = default;
        DiagonalOperator& operator=( DiagonalOperator&& vec ) = default; 

        explicit DiagonalOperator( int, Float s );
        explicit DiagonalOperator( int, const FloatVector& dia );
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
        virtual void print( std::ostream& ) const override;

        FloatVector& getdiagonal();
        const FloatVector& getdiagonal() const;
        
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
    
    return DiagonalOperator( left.getdimout(), 
                             FloatVector( leftdia.getdimension(), 
                                          [&](int d) -> Float { 
                                            assert( 0 <= d && d < dimension ); 
                                            return leftdia[d] * rightdia[d];
                                          }
                                        )
                           );
}  
  
  
  
  

#endif
