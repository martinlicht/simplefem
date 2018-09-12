

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Index Ranges" << std::endl;
        
        if( true ) {
          
          cout << "Test empty index ranges" << std::endl;
          
          IndexRange irE( 3, 2 );
          
          assert( irE.isempty() );
          assert( irE.cardinality() == 0 );
          assert( irE.cardinality() );
          
          for( int i : irE )
            cout << "Output :" << i << endl;
          
        
        } else {
          
          cout << "Test empty index range" << std::endl;
          
        } 
          
        if( true ) {
          
          cout << "Test non-empty index ranges" << std::endl;
          
          IndexRange irA( 3, 7 );
          IndexRange irB( 5, 5 );
          
          assert( irA.cardinality() == 5 );
          
          assert( irA.cardinality() == irA.cardinality() );
          
          assert( !irA.contains(2) );
          assert( irA.contains(3) );
          assert( irA.contains(4) );
          assert( irA.contains(5) );
          assert( irA.contains(6) );
          assert( irA.contains(7) );
          assert( !irA.contains(8) );
          
          assert( irB.cardinality() == irB.cardinality() );
          
          assert( irB.cardinality() == 1 );
          
          assert( irB.contains(5) );
          assert( !irB.contains(4) );
          assert( !irB.contains(6) );
          
          cout << "Test indexing in non-empty index ranges" << std::endl;
          
          assert( irB.element2position(5) == 0 );
          assert( irB.position2element(0) == 5 );
          assert( irA.element2position(3) == 0 );
          assert( irA.element2position(5) == 2 );
          assert( irA.element2position(7) == 4 );
          assert( irA.position2element(0) == 3 );
          assert( irA.position2element(2) == 5 );
          assert( irA.position2element(4) == 7 );
          
          
        }
        
        if( true ) {
          
          cout << "Test Index Range iterators" << std::endl;
          
          IndexRange irC( 2, 5 );
          
          assert( !irC.isempty() );
          
          cout << "For each loop " << std::endl;
          
          for( int i : irC )
            cout << "Output :" << i << endl;
            
          cout << "Classical For loop " << std::endl;
          
          for( IndexRange::IndexRangeIterator iri = irC.begin(); iri != irC.end(); iri++ )
            cout << "Output :" << *iri << endl;
            
          cout << "While Loop " << std::endl;
          
          IndexRange::IndexRangeIterator iri = irC.begin();
          while( iri != irC.end() )
            cout << "Output :" << *(iri++) << endl;
            
          
          
        } else {
          
          cout << "Test empty index range" << std::endl;
          
        } 
          
        if( false ) {
        
          cout << "Test invalid index ranges" << std::endl;
        
          IndexRange dummy( 0, std::numeric_limits<int>::max() );
          cout << dummy << " " << std::numeric_limits<int>::max();
          cout << "FAIL: Invalid range created" << std::endl;
        
        } else {
          
          cout << "No Test for invalid index ranges" << std::endl;
          
        }
          
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
