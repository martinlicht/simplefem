#ifndef INCLUDEGUARD_OPERATOR_COMPOSEDOPERATORS_HPP
#define INCLUDEGUARD_OPERATOR_COMPOSEDOPERATORS_HPP

#include <cassert>
#include <memory>
#include <ostream>
#include <utility>

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "simpleoperators.hpp"


/************************
****
****  Classes for Composition of Operators TODO
****  - instantiates LinearOperator
****  
************************/


class ProxyOperator final
: public LinearOperator 
{

    public:

        explicit ProxyOperator() = delete;
        explicit ProxyOperator( const ProxyOperator& ) = delete;
        explicit ProxyOperator( ProxyOperator&& ) = delete;
        ProxyOperator& operator=( const ProxyOperator& vec ) = delete;
        ProxyOperator& operator=( ProxyOperator&& vec ) = delete; 

        
        explicit ProxyOperator( const LinearOperator& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), op( op ) { 
//             LOG << "Proxy created" << ""; 
        }
        
        virtual ~ProxyOperator() { 
//             LOG << "Proxy destroyed" << "";
        }

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<ProxyOperator> cloned = std::make_shared<ProxyOperator>( op );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<ProxyOperator> heir = std::make_unique<ProxyOperator>( op );
            return heir;
        }

        virtual void check() const override { 
            op.check();
            assert( getdimout() == op.getdimout() );
            assert( getdimin()  == op.getdimin()  );
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Proxy Operator" << std::endl;
        }
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
            check();
            src.check();
            dest.check();
            
            assert( getdimin() == src.getdimension() );
            assert( getdimout() == dest.getdimension() );
            
//             LOG << "call";
//             op.print( std::cout );
            op.apply( dest, src, scaling );    
//             LOG << "done";
        }

    private:

        const LinearOperator& op;
    
};




class ComposedOperator
: public LinearOperator 
{

    public:
    
        explicit ComposedOperator() = delete;
        explicit ComposedOperator( const ComposedOperator& ) = delete;
        explicit ComposedOperator( ComposedOperator&& ) = delete;
        ComposedOperator& operator=( const ComposedOperator& vec ) = delete;
        ComposedOperator& operator=( ComposedOperator&& vec ) = delete; 

        
        explicit ComposedOperator( int dimout, int dimin, std::unique_ptr<LinearOperator>&& pl, std::unique_ptr<LinearOperator>&& pr )
        : LinearOperator( dimout, dimin ), left( std::move(pl) ), right( std::move(pr) ) 
        {
            ComposedOperator::check();
        } 
        
    public:

        explicit ComposedOperator( int dimout, int dimin, const LinearOperator& L, const LinearOperator& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
//             assert( L.getdimin() == R.getdimout() );
            left  = std::make_unique<ProxyOperator>(L); 
            right = std::make_unique<ProxyOperator>(R);
            ComposedOperator::check();
//             LOG << "Composed Operator created LL" << ""; 
        }
        
        explicit ComposedOperator( int dimout, int dimin, const LinearOperator& L, LinearOperator&& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
//             assert( L.getdimin() == R.getdimout() );
            left  = std::make_unique<ProxyOperator>(L); 
            right = std::move(R).get_unique_pointer_to_heir();
            ComposedOperator::check();
//             LOG << "Composed Operator created LR" << ""; 
        }
        
        explicit ComposedOperator( int dimout, int dimin, LinearOperator&& L, const LinearOperator& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
//             assert( L.getdimin() == R.getdimout() );
            left  = std::move(L).get_unique_pointer_to_heir();
            right = std::make_unique<ProxyOperator>(R);
            ComposedOperator::check();
//             LOG << "Composed Operator created RL" << ""; 
        }
        
        explicit ComposedOperator( int dimout, int dimin, LinearOperator&& L, LinearOperator&& R )
        : LinearOperator( dimout, dimin ), left( nullptr ), right( nullptr ) 
        {
//             assert( L.getdimin() == R.getdimout() );
            left  = std::move(L).get_unique_pointer_to_heir();
            right = std::move(R).get_unique_pointer_to_heir();
            ComposedOperator::check();
//             LOG << "Composed Operator created RR" << ""; 
        }
        
        virtual ~ComposedOperator() { 
            //ComposedOperator::check(); // explicitly disabled, might be in moved-from state 
//             LOG << "Composed Operator destroyed" << "";
        }

        virtual void check() const override { 
            LinearOperator::check(); 
            assert( left != nullptr );
            assert( right != nullptr );
            left->check(); right->check(); 
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








class ProduktOperator final
: public ComposedOperator 
{

    public:
    
        explicit ProduktOperator() = delete;
        explicit ProduktOperator( const ProduktOperator& ) = delete;
        explicit ProduktOperator( ProduktOperator&& ) = delete;
        ProduktOperator& operator=( const ProduktOperator& vec ) = delete;
        ProduktOperator& operator=( ProduktOperator&& vec ) = delete;

        explicit ProduktOperator( int dimout, int dimin, std::unique_ptr<LinearOperator>&& pl, std::unique_ptr<LinearOperator>&& pr )
        : ComposedOperator( dimout, dimin, std::move(pl), std::move(pr) ) 
        {
            ProduktOperator::check();
        } 
        
    public:

        explicit ProduktOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            ProduktOperator::check();
        }
        
        explicit ProduktOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            ProduktOperator::check();
        }
        
        explicit ProduktOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            ProduktOperator::check();
        }
        
        explicit ProduktOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            ProduktOperator::check();
        }
        
        virtual ~ProduktOperator()
        {
//             ProduktOperator::check();
        }

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            check();
            auto lp = std::move( *(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto rp = std::move( *(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            std::shared_ptr<ProduktOperator> clone = std::make_shared<ProduktOperator>( getdimout(), getdimin(), std::move(lp), std::move(rp) );
            clone->check();
            return clone;
//             LOG << "Composed Operator created (clone)" << ""; 
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            check();
            std::unique_ptr<ProduktOperator> heir = std::make_unique<ProduktOperator>( getdimout(), getdimin(), std::move(left), std::move(right) );
            heir->check();
//             LOG << "Composed Operator created (heir)" << ""; 
            return heir;
        }
        
        virtual void check() const override { 
            ComposedOperator::check();
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimin() == right->getdimout() );
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Produkt Operator" << nl;
            ComposedOperator::print( os );
        }
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override {
            dest = scaling * left->apply( right->apply(add) );
        }

};

inline ProduktOperator operator*( const LinearOperator& left, const LinearOperator& right )
{ return ProduktOperator( left, right ); }

inline ProduktOperator operator*( const LinearOperator& left, LinearOperator&& right )
{ return ProduktOperator( left, std::move(right) ); }

inline ProduktOperator operator*( LinearOperator&& left, const LinearOperator& right )
{ return ProduktOperator( std::move(left), right ); }

inline ProduktOperator operator*( LinearOperator&& left, LinearOperator&& right )
{ return ProduktOperator( std::move(left), std::move(right) ); }
  








class SummOperator final
: public ComposedOperator 
{

    public:
    
        explicit SummOperator() = delete;
        explicit SummOperator( const SummOperator& ) = delete;
        explicit SummOperator( SummOperator&& ) = delete;
        SummOperator& operator=( const SummOperator& vec ) = delete;
        SummOperator& operator=( SummOperator&& vec ) = delete; 

        explicit SummOperator( int dimout, int dimin, std::unique_ptr<LinearOperator>&& pl, std::unique_ptr<LinearOperator>&& pr )
        : ComposedOperator( dimout, dimin, std::move(pl), std::move(pr) ) 
        {
            SummOperator::check();
        } 
        
    public:

        explicit SummOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            SummOperator::check();
        }
        
        explicit SummOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            SummOperator::check();
        }
        
        explicit SummOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            SummOperator::check();
        }
        
        explicit SummOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            SummOperator::check();
        }
        
        virtual ~SummOperator()
        {
//             SummOperator::check();
        }
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            check();
            auto lp = std::move( *(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto rp = std::move( *(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            std::shared_ptr<SummOperator> clone = std::make_shared<SummOperator>( getdimout(), getdimin(), std::move(lp), std::move(rp) );
            clone->check();
            return clone;
//             LOG << "Composed Operator created (clone)" << ""; 
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            check();
            std::unique_ptr<SummOperator> heir = std::make_unique<SummOperator>( getdimout(), getdimin(), std::move(left), std::move(right) );
            heir->check();
//             LOG << "Composed Operator created (heir)" << ""; 
            return heir;
        }
        
        virtual void check() const override
        { 
            ComposedOperator::check();
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
        }
        
        virtual void print( std::ostream& os ) const override
        { 
            os << "Print Summ Operator" << nl;
            ComposedOperator::print( os );
        }
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) + scaling * right->apply( add );
        }

};

inline SummOperator operator+( const LinearOperator& left, const LinearOperator& right )
{ return SummOperator( left, right ); }

inline SummOperator operator+( const LinearOperator& left, LinearOperator&& right )
{ return SummOperator( left, std::move(right) ); }

inline SummOperator operator+( LinearOperator&& left, const LinearOperator& right )
{ return SummOperator( std::move(left), right ); }

inline SummOperator operator+( LinearOperator&& left, LinearOperator&& right )
{ return SummOperator( std::move(left), std::move(right) ); }
 






class DiffOperator final
: public ComposedOperator 
{

    public:
    
        explicit DiffOperator() = delete;
        explicit DiffOperator( const DiffOperator& ) = delete;
        explicit DiffOperator( DiffOperator&& ) = delete;
        DiffOperator& operator=( const DiffOperator& vec ) = delete;
        DiffOperator& operator=( DiffOperator&& vec ) = delete; 

        explicit DiffOperator( int dimout, int dimin, std::unique_ptr<LinearOperator>&& pl, std::unique_ptr<LinearOperator>&& pr )
        : ComposedOperator( dimout, dimin, std::move(pl), std::move(pr) ) 
        {
            DiffOperator::check();
        } 
        
    public:

        explicit DiffOperator( const LinearOperator& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, R ) 
        {
            DiffOperator::check();
        }
        
        explicit DiffOperator( const LinearOperator& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), L, std::move(R) ) 
        {
            DiffOperator::check();
        }
        
        explicit DiffOperator( LinearOperator&& L, const LinearOperator& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), R ) 
        {
            DiffOperator::check();
        }
        
        explicit DiffOperator( LinearOperator&& L, LinearOperator&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::move(L), std::move(R) ) 
        {
            DiffOperator::check();
        }
        
        virtual ~DiffOperator()
        {
//             DiffOperator::check();
        }
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            check();
            auto lp = std::move( *(left ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto rp = std::move( *(right->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            std::shared_ptr<DiffOperator> clone = std::make_shared<DiffOperator>( getdimout(), getdimin(), std::move(lp), std::move(rp) );
            clone->check();
            return clone;
//             LOG << "Composed Operator created (clone)" << ""; 
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            check();
            std::unique_ptr<DiffOperator> heir = std::make_unique<DiffOperator>( getdimout(), getdimin(), std::move(left), std::move(right) );
            heir->check();
//             LOG << "Composed Operator created (heir)" << ""; 
            return heir;
        }
        
        virtual void check() const override
        { 
            ComposedOperator::check();
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
        }
        
        virtual void print( std::ostream& os ) const override
        { 
            os << "Print Diff Operator" << nl;
            ComposedOperator::print( os );
        }
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) - scaling * right->apply( add );
        }

};

inline DiffOperator operator-( const LinearOperator& left, const LinearOperator& right )
{ return DiffOperator( left, right ); }

inline DiffOperator operator-( const LinearOperator& left, LinearOperator&& right )
{ return DiffOperator( left, std::move(right) ); }

inline DiffOperator operator-( LinearOperator&& left, const LinearOperator& right )
{ return DiffOperator( std::move(left), right ); }

inline DiffOperator operator-( LinearOperator&& left, LinearOperator&& right )
{ return DiffOperator( std::move(left), std::move(right) ); }
 

 
















inline ProduktOperator operator-( LinearOperator&& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), -1. ), std::move(op) );
}

inline ProduktOperator operator-( LinearOperator& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), -1. ), op );
}

inline ProduktOperator operator+( LinearOperator&& op )
{
    return ProduktOperator( IdentityOperator( op.getdimout() ), std::move(op) );
}

inline ProduktOperator operator+( LinearOperator& op )
{
    return ProduktOperator( IdentityOperator( op.getdimout() ), op );
}

inline ProduktOperator operator*( Float s, LinearOperator&& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), s ), std::move(op) );
}

inline ProduktOperator operator*( Float s, LinearOperator& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), s ), op );
}













class Block2x2Operator
: public LinearOperator 
{

    public:
    
        explicit Block2x2Operator() = delete;
        explicit Block2x2Operator( const Block2x2Operator& ) = delete;
        explicit Block2x2Operator( Block2x2Operator&& ) = delete;
        Block2x2Operator& operator=( const Block2x2Operator& vec ) = delete;
        Block2x2Operator& operator=( Block2x2Operator&& vec ) = delete; 

        
        explicit Block2x2Operator( int dimout, int dimin,
                                   std::unique_ptr<LinearOperator>&& pul, std::unique_ptr<LinearOperator>&& pur, 
                                   std::unique_ptr<LinearOperator>&& pll, std::unique_ptr<LinearOperator>&& plr 
                                 )
        : LinearOperator( dimout, dimin ), 
        upperleft( std::move(pul) ), upperright( std::move(pur) ),
        lowerleft( std::move(pll) ), lowerright( std::move(plr) ) 
        {
            Block2x2Operator::check();
        } 
        
    public:

//         left  = std::make_unique<ProxyOperator>(L); 
//         right = std::move(R).get_unique_pointer_to_heir();
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL, const LinearOperator&  UR, 
                                   const LinearOperator&  LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL, const LinearOperator&  UR, 
                                   const LinearOperator&  LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL, const LinearOperator&  UR, 
                                         LinearOperator&& LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL, const LinearOperator&  UR, 
                                         LinearOperator&& LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL, const LinearOperator&  UR, 
                                   const LinearOperator&  LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL, const LinearOperator&  UR, 
                                   const LinearOperator&  LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL, const LinearOperator&  UR, 
                                         LinearOperator&& LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL, const LinearOperator&  UR, 
                                         LinearOperator&& LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::make_unique<ProxyOperator>(UR);
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL,       LinearOperator&& UR, 
                                   const LinearOperator&  LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL,       LinearOperator&& UR, 
                                   const LinearOperator&  LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL,       LinearOperator&& UR, 
                                         LinearOperator&& LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                   const LinearOperator&  UL,       LinearOperator&& UR, 
                                         LinearOperator&& LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::make_unique<ProxyOperator>(UL); 
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL,       LinearOperator&& UR, 
                                   const LinearOperator&  LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL,       LinearOperator&& UR, 
                                   const LinearOperator&  LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::make_unique<ProxyOperator>(LL); 
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL,       LinearOperator&& UR, 
                                         LinearOperator&& LL, const LinearOperator&  LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::make_unique<ProxyOperator>(LR);
            Block2x2Operator::check();
        }
        
        explicit Block2x2Operator( int dimout, int dimin, 
                                         LinearOperator&& UL,       LinearOperator&& UR, 
                                         LinearOperator&& LL,       LinearOperator&& LR 
                                 )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft( nullptr ), upperright( nullptr ),
        lowerleft( nullptr ), lowerright( nullptr ) 
        {
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            upperleft  = std::move(UL).get_unique_pointer_to_heir();
            upperright = std::move(UR).get_unique_pointer_to_heir();
            lowerleft  = std::move(LL).get_unique_pointer_to_heir();
            lowerright = std::move(LR).get_unique_pointer_to_heir();
            Block2x2Operator::check();
        }
        
        
        
        
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            check();
            auto ulp = std::move( *(upperleft ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto urp = std::move( *(upperright->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto llp = std::move( *(lowerleft ->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            auto lrp = std::move( *(lowerright->get_shared_pointer_to_clone())).get_unique_pointer_to_heir();
            std::shared_ptr<Block2x2Operator> clone = std::make_shared<Block2x2Operator>( getdimout(), getdimin(), std::move(ulp), std::move(urp), std::move(llp), std::move(lrp) );
            clone->check();
            return clone;
//             LOG << "Composed Operator created (clone)" << ""; 
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            check();
            std::unique_ptr<Block2x2Operator> heir = std::make_unique<Block2x2Operator>( getdimout(), getdimin(), std::move(upperleft), std::move(upperright), std::move(lowerleft), std::move(lowerright) );
            heir->check();
//             LOG << "Composed Operator created (heir)" << ""; 
            return heir;
        }

        
        
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override {
            
            assert( dest.getdimension() == upperleft->getdimout()  + lowerleft->getdimout()  );
            assert( dest.getdimension() == upperright->getdimout() + lowerright->getdimout() );
            assert(  add.getdimension() == upperleft->getdimin()   + upperright->getdimin()  );
            assert(  add.getdimension() == lowerleft->getdimin()   + lowerright->getdimin()  );
            
            
            const auto left  = add.getslice( 0,                     upperleft ->getdimin() );
            const auto right = add.getslice( upperleft->getdimin(), upperright->getdimin() );
            
            const auto upper = (*upperleft) * left + (*upperright) * right;
            const auto lower = (*lowerleft) * left + (*lowerright) * right;
            
            dest.setslice( 0,                    upper );
            dest.setslice( upper.getdimension(), lower );
            
            dest.scale( scaling );
        }

        
        
        
        virtual ~Block2x2Operator() {}

        virtual void check() const override { 
            LinearOperator::check(); 
            assert( upperleft  != nullptr );
            assert( upperright != nullptr );
            assert( lowerleft  != nullptr );
            assert( lowerright != nullptr );
            upperleft->check();
            upperright->check(); 
            lowerleft->check();
            lowerright->check(); 
        }
        
        virtual void print( std::ostream& os ) const override { 
            os << "Print Composed Operator" << nl;
            upperleft->print(os); 
            upperright->print(os); 
            lowerleft->print(os); 
            lowerright->print(os); 
        }
        
    private:

        std::unique_ptr<LinearOperator> upperleft;
        std::unique_ptr<LinearOperator> upperright;
        std::unique_ptr<LinearOperator> lowerleft;
        std::unique_ptr<LinearOperator> lowerright;
        
};



#endif
