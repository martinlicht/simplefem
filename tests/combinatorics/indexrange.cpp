

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"



using namespace std;

int main()
{
    LOG << "Unit Test for Index Ranges" << nl;
    
    if( true ) 
    {
        
        LOG << "Test empty index ranges" << nl;
        
        IndexRange irE1(  3,  2 );
        IndexRange irE2(  5,  1 );
        IndexRange irE3( -2, -3 );
        
        assert( irE1.isempty() );
        assert( irE2.isempty() );
        assert( irE3.isempty() );
        
        assert( irE1.cardinality() == 0 );
        assert( irE2.cardinality() == 0 );
        assert( irE3.cardinality() == 0 );
        
        for( int i : irE1 ) assert( false );
        for( int i : irE2 ) assert( false );
        for( int i : irE3 ) assert( false );
        
        assert( irE1 == irE2 );
        assert( irE1 == irE3 );
        assert( irE2 == irE3 );
        
    } 
        
    if( true ) 
    {
        
        LOG << "Test non-empty index ranges" << nl;
        
        IndexRange irA( 3, 7 );
        IndexRange irB( 5, 5 );
        IndexRange irC(-2,-2 );
        IndexRange irD( 4, 9 );
        
        assert( not irA.isempty() );
        assert( not irB.isempty() );
        assert( not irC.isempty() );
        assert( not irD.isempty() );
        
        assert( irA.cardinality() == 5 );
        assert( irB.cardinality() == 1 );
        assert( irC.cardinality() == 1 );
        assert( irD.cardinality() == 6 );
        
        assert( irA != irB );
        assert( irA != irC );
        assert( irA != irD );
        assert( irB != irC );
        assert( irB != irD );
        assert( irC != irD );
        
        assert( !irA.contains(2) );
        assert( irA.contains(3) );
        assert( irA.contains(4) );
        assert( irA.contains(5) );
        assert( irA.contains(6) );
        assert( irA.contains(7) );
        assert( !irA.contains(8) );
        
        assert( irB.contains(5) );
        assert( !irB.contains(4) );
        assert( !irB.contains(6) );
        
        LOG << "Test indexing in non-empty index ranges" << nl;
        
        assert( irB.element2position(5) == 0 );
        assert( irB.position2element(0) == 5 );
        assert( irA.element2position(3) == 0 );
        assert( irA.element2position(5) == 2 );
        assert( irA.element2position(7) == 4 );
        assert( irA.position2element(0) == 3 );
        assert( irA.position2element(2) == 5 );
        assert( irA.position2element(4) == 7 );
        
    }
    
    if( true )
    {
        
        LOG << "Test Index Range iterators" << nl;
        
        IndexRange irC( 2, 5 );
        
        assert( not irC.isempty() );
        
        LOG << "For each loop " << nl;
        for( int i : irC )
            LOG << i << space;
        LOG << nl;

        LOG << "Classical For loop " << nl;
        for( IndexRange::ConstIterator iri = irC.begin(); iri != irC.end(); iri++ )
            LOG << *iri << space;
        LOG << nl;
        
        LOG << "While Loop " << nl;
        IndexRange::ConstIterator iri = irC.begin();
        while( iri != irC.end() )
            LOG << *(iri++) << space << nl;
        LOG << nl;
        
    }
        
    if( false )
    {
    
        LOG << "Test invalid index ranges" << nl;
    
        IndexRange dummy( 0, std::numeric_limits<int>::max() );
        LOG << dummy << " " << std::numeric_limits<int>::max();
        LOG << "FAIL: Invalid range created" << nl;
    
    }
        
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
