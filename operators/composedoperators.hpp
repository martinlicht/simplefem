#ifndef INCLUDEGUARD_OPERATOR_COMPOSED
#define INCLUDEGUARD_OPERATOR_COMPOSED

#include <memory>

#include "../basic.hpp"
#include "linearoperator.hpp"


/************************
****
****  Class for Scalings 
****  - instantiates LinearOperator
****  
************************/


class ProxyOperator:
public LinearOperator /* every scaling operation is a linear operator */
{

    public:

        explicit ProxyOperator( const LinearOperator& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), op( op ) { 
            std::cout << "Proxy created" << nl; 
        };
        
        virtual ~ProxyOperator() { 
            std::cout << "Proxy destroyed" << nl;
        };

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<ProxyOperator> clone = std::make_shared<ProxyOperator>( op );
            return clone;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<ProxyOperator> heir = std::make_unique<ProxyOperator>( op );
            return heir;
        }

        virtual void check() const override { 
            op.check();
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Proxy Operator" << std::endl;
        }
        
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
            check();
            src.check();
            dest.check();
            
            assert( getdimin() == src.getdimension() );
            assert( getdimout() == dest.getdimension() );
            
            op.apply( dest, src, scaling );    
        }

    private:

        const LinearOperator& op;
    
};




class ComposedOperator:
public LinearOperator 
{

    public:
    
        explicit ComposedOperator( int dimout, int dimin )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {}; 
        
    public:

        explicit ComposedOperator( int dimout, int dimin, const LinearOperator& L, const LinearOperator& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
            assert( L.getdimin() == R.getdimout() );
            left  = std::make_unique<ProxyOperator>(L); 
            right = std::make_unique<ProxyOperator>(R);
            check();
            std::cout << "Composed Operator created LL" << nl; 
        };
        
        explicit ComposedOperator( int dimout, int dimin, const LinearOperator& L, LinearOperator&& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
            assert( L.getdimin() == R.getdimout() );
            left  = std::make_unique<ProxyOperator>(L); 
            right = std::move(R).get_unique_pointer_to_heir();
            check();
            std::cout << "Composed Operator created LR" << nl; 
        };
        
        explicit ComposedOperator( int dimout, int dimin, LinearOperator&& L, const LinearOperator& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
            assert( L.getdimin() == R.getdimout() );
            left  = std::move(L).get_unique_pointer_to_heir();
            right = std::make_unique<ProxyOperator>(R);
            check();
            std::cout << "Composed Operator created RL" << nl; 
        };
        
        explicit ComposedOperator( int dimout, int dimin, LinearOperator&& L, LinearOperator&& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
            assert( L.getdimin() == R.getdimout() );
            left  = std::move(L).get_unique_pointer_to_heir();
            right = std::move(R).get_unique_pointer_to_heir();
            check();
            std::cout << "Composed Operator created RR" << nl; 
        };
        
        virtual ~ComposedOperator() { 
            std::cout << "Composed Operator destroyed" << nl;
        };

        virtual void check() const override { 
            LinearOperator::check(); left->check(); right->check(); 
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Composed Operator" << nl;
            left->print(os); 
            right->print(os); 
        }
        
    protected:

        std::unique_ptr<LinearOperator> left;
        std::unique_ptr<LinearOperator> right;
        
};








class ProduktOperator:
public ComposedOperator 
{

    public:
    
        explicit ProduktOperator( int dimout, int dimin )
        : ComposedOperator( dimout, dimin ) 
        {}; 
        
    public:

        explicit ProduktOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            check();
        };
        
        explicit ProduktOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            check();
        };
        
        explicit ProduktOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            check();
        };
        
        explicit ProduktOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            check();
        };
        
        virtual ~ProduktOperator() = default;
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<ProduktOperator> clone = std::make_shared<ProduktOperator>( getdimout(), getdimin() );
            clone->left  = std::move(*(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->right = std::move(*(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->check();
            return clone;
            std::cout << "Composed Operator created (clone)" << nl; 
        };
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<ProduktOperator> heir = std::make_unique<ProduktOperator>( getdimout(), getdimin() );
            heir->left  = std::move(left);
            heir->right = std::move(right);
            heir->check();
            std::cout << "Composed Operator created (heir)" << nl; 
            return heir;
        };
        
        virtual void check() const override { 
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimin() == right->getdimout() );
            ComposedOperator::check();
        };
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Produkt Operator" << nl;
            ComposedOperator::print( os );
        };
        
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override {
            dest = scaling * left->apply( right->apply(add) );
        };

};

inline ProduktOperator operator*( const LinearOperator& left, const LinearOperator& right )
{ return ProduktOperator( left, right ); }

inline ProduktOperator operator*( const LinearOperator& left, LinearOperator&& right )
{ return ProduktOperator( left, std::move(right) ); }

inline ProduktOperator operator*( LinearOperator&& left, const LinearOperator& right )
{ return ProduktOperator( std::move(left), right ); }

inline ProduktOperator operator*( LinearOperator&& left, LinearOperator&& right )
{ return ProduktOperator( std::move(left), std::move(right) ); }
  








class SummOperator:
public ComposedOperator 
{

    public:
    
        explicit SummOperator( int dimout, int dimin )
        : ComposedOperator( dimout, dimin ) 
        {}; 
        
    public:

        explicit SummOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            check();
        };
        
        explicit SummOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            check();
        };
        
        explicit SummOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            check();
        };
        
        explicit SummOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            check();
        };
        
        virtual ~SummOperator() = default;
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<SummOperator> clone = std::make_shared<SummOperator>( getdimout(), getdimin() );
            clone->left  = std::move(*(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->right = std::move(*(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->check();
            return clone;
            std::cout << "Composed Operator created (clone)" << nl; 
        };
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<SummOperator> heir = std::make_unique<SummOperator>( getdimout(), getdimin() );
            heir->left  = std::move(left);
            heir->right = std::move(right);
            heir->check();
            std::cout << "Composed Operator created (heir)" << nl; 
            return heir;
        };
        
        virtual void check() const override
        { 
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
            ComposedOperator::check();
        };
        
        virtual void print( std::ostream& os ) const override
        { 
            os << "Print Summ Operator" << nl;
            ComposedOperator::print( os );
        };
        
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) + scaling * right->apply( add );
        };

};

inline SummOperator operator+( const LinearOperator& left, const LinearOperator& right )
{ return SummOperator( left, right ); }

inline SummOperator operator+( const LinearOperator& left, LinearOperator&& right )
{ return SummOperator( left, std::move(right) ); }

inline SummOperator operator+( LinearOperator&& left, const LinearOperator& right )
{ return SummOperator( std::move(left), right ); }

inline SummOperator operator+( LinearOperator&& left, LinearOperator&& right )
{ return SummOperator( std::move(left), std::move(right) ); }
 






class DiffOperator:
public ComposedOperator 
{

    public:
    
        explicit DiffOperator( int dimout, int dimin )
        : ComposedOperator( dimout, dimin ) 
        {}; 
        
    public:

        explicit DiffOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            check();
        };
        
        explicit DiffOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            check();
        };
        
        explicit DiffOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            check();
        };
        
        explicit DiffOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            check();
        };
        
        virtual ~DiffOperator() = default;
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<DiffOperator> clone = std::make_shared<DiffOperator>( getdimout(), getdimin() );
            clone->left  = std::move(*(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->right = std::move(*(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            clone->check();
            return clone;
            std::cout << "Composed Operator created (clone)" << nl; 
        };
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<DiffOperator> heir = std::make_unique<DiffOperator>( getdimout(), getdimin() );
            heir->left  = std::move(left);
            heir->right = std::move(right);
            heir->check();
            std::cout << "Composed Operator created (heir)" << nl; 
            return heir;
        };
        
        virtual void check() const override
        { 
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
            ComposedOperator::check();
        };
        
        virtual void print( std::ostream& os ) const override
        { 
            os << "Print Diff Operator" << nl;
            ComposedOperator::print( os );
        };
        
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) - scaling * right->apply( add );
        };

};

inline DiffOperator operator-( const LinearOperator& left, const LinearOperator& right )
{ return DiffOperator( left, right ); }

inline DiffOperator operator-( const LinearOperator& left, LinearOperator&& right )
{ return DiffOperator( left, std::move(right) ); }

inline DiffOperator operator-( LinearOperator&& left, const LinearOperator& right )
{ return DiffOperator( std::move(left), right ); }

inline DiffOperator operator-( LinearOperator&& left, LinearOperator&& right )
{ return DiffOperator( std::move(left), std::move(right) ); }
 

 

  

#endif
