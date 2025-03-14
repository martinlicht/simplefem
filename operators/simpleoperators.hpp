#ifndef INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS_HPP
#define INCLUDEGUARD_OPERATOR_SIMPLEOPERATORS_HPP

#include "../base/include.hpp"
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

        /* Constructors */
        
        explicit IdentityOperator( int n );
        
        /* standard interface */
        
        IdentityOperator()                                        = delete;
        IdentityOperator( const IdentityOperator& )               = default;
        IdentityOperator( IdentityOperator&& )                    noexcept = default;
        IdentityOperator& operator=( const IdentityOperator& op ) = default;
        IdentityOperator& operator=( IdentityOperator&& op )      noexcept = default; 
        virtual ~IdentityOperator()                               noexcept = default;

        /* standard methods for operators */
        
        virtual void check() const override;
        
        virtual std::string text() const override;

        /* OTHER METHODS */
        
        virtual IdentityOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
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
****  Class for Zero Operators 
****  
************************/


class ZeroOperator final
: public LinearOperator
{

    public:

        /* Constructors */
        
        explicit ZeroOperator( int dim );
        explicit ZeroOperator( int dimout, int dimin );
        
        /* standard interface */
        
        ZeroOperator()                                    = delete;
        ZeroOperator( const ZeroOperator& )               = default;
        ZeroOperator( ZeroOperator&& )                    noexcept = default;
        ZeroOperator& operator=( const ZeroOperator& op ) = default;
        ZeroOperator& operator=( ZeroOperator&& op )      noexcept = default; 
        virtual ~ZeroOperator()                           noexcept = default;

        /* standard methods for operators */
        
        virtual void check() const override;
        
        virtual std::string text() const override;

        /* OTHER METHODS */
        
        virtual ZeroOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
    
};
  
  







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

        /* Constructors */
        
        explicit ScalingOperator( int, Float s );
        /* standard methods for operators */
        ScalingOperator()                                       = delete;
        ScalingOperator( const ScalingOperator& )               = default;
        ScalingOperator( ScalingOperator&& )                    noexcept = default;
        ScalingOperator& operator=( const ScalingOperator& op ) = default;
        ScalingOperator& operator=( ScalingOperator&& op )      noexcept = default; 
        virtual ~ScalingOperator()                              noexcept = default;

        
        /* standard interface */
        
        virtual void check() const override;
        
        virtual std::string text() const override;

        /* OTHER METHODS */
        
        virtual ScalingOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        Float getscaling() const;
        void setscaling( Float s );

        using LinearOperator::apply; // import any 'apply' into the derived class' methods
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

        /* Constructors */
        
        explicit DiagonalOperator( int, Float s );
        explicit DiagonalOperator( const FloatVector& dia );
        explicit DiagonalOperator( FloatVector&& dia );
        explicit DiagonalOperator( int, const ScalingOperator& scaling );
        explicit DiagonalOperator( int, const std::function<Float(int)>& );
        
        /* standard methods for operators */
        
        DiagonalOperator()                                        = delete;
        DiagonalOperator( const DiagonalOperator& )               = default;
        DiagonalOperator( DiagonalOperator&& )                    noexcept = default;
        DiagonalOperator& operator=( const DiagonalOperator& op ) = default;
        DiagonalOperator& operator=( DiagonalOperator&& op )      noexcept = default; 
        virtual ~DiagonalOperator()                               noexcept = default;
        
        /* standard interface */
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual std::string text( const bool embellish ) const;

        /* OTHER METHODS */
        
        virtual DiagonalOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        

        FloatVector& getDiagonal();
        const FloatVector& getDiagonal() const;
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
        
        DiagonalOperator sqrt() const;

    private:

        FloatVector diagonal;
            
};
  
  

inline DiagonalOperator operator*( const DiagonalOperator& left, const DiagonalOperator& right )
{
    left.check();
    right.check();
    
    assert( left.getdimin() == right.getdimout() );
    assert( left.getdimout() == right.getdimin() );
    
    const FloatVector& leftdia = left.getDiagonal();
    const FloatVector& rightdia = right.getDiagonal();
    // const int dimension = leftdia.getdimension();
    
    assert( leftdia.getdimension() == rightdia.getdimension() );
    
    return DiagonalOperator( FloatVector( leftdia.getdimension(), 
                                          [&](int d) -> Float { 
                                            assert( 0 <= d && d < leftdia.getdimension() ); 
                                            return leftdia[d] * rightdia[d];
                                          }
                                        )
                           );
}  
  
  
  
  









class LambdaOperator final
: public LinearOperator
{

    public:

        /* Constructors */
        
        explicit LambdaOperator( int dim, const std::function<FloatVector(const FloatVector&)>& );
        explicit LambdaOperator( int dimout, int dimin, const std::function<FloatVector(const FloatVector&)>& );
        
        /* standard interface */
        
        LambdaOperator()                                      = delete;
        LambdaOperator( const LambdaOperator& )               = default;
        LambdaOperator( LambdaOperator&& )                    noexcept = default;
        LambdaOperator& operator=( const LambdaOperator& op ) = delete;
        LambdaOperator& operator=( LambdaOperator&& op )      noexcept = delete; 
        virtual ~LambdaOperator()                             noexcept = default;

        /* standard methods for operators */
        
        virtual void check() const override;
        
        virtual std::string text() const override;

        /* OTHER METHODS */
        
        virtual LambdaOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;

    private:

        const std::function<FloatVector(const FloatVector&)> func;
    
};


















#endif
