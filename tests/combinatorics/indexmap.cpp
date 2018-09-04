

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Index Mapping" << endl;

    const IndexRange irA( 2, 5 );
    const IndexRange irB( 3, 5 );
    const IndexRange irC( -2, 2 );
    const IndexRange irD( 0, 7 );

    {
      
      cout << "Test Identity of usual interval" << endl;
        
      const IndexMap id  = identityIndexMap( irA );
      cout << id << endl;
      
      assert( id == identityIndexMap( irA.min(), irA.max() ) );
      
      for( IndexRange::IndexRangeIterator iri : irA ) {
        assert( id.rangecontains( *iri ) );
        assert( id.rangeposition( *iri ) == *iri );
      }
      
      assert( id.isbijective()  );
      assert( id.isinjective()  );
      assert( id.issurjective() );
      assert( id.isstrictlyascending() );
      
      assert( id.comparablewith( id ) );
      assert( id.equals( id ) );
      assert( ! id.less( id ) );
      
    }
    
    {
      
      cout << "Test Empty Index Map" << endl;
        
      const IndexMap leer  = identityIndexMap( IndexRange( 1, 0 ) );
      cout << leer << endl;
      
      assert( leer == identityIndexMap( 1, 0 ) );
      
      assert( leer.isbijective()  );
      assert( leer.isinjective()  );
      assert( leer.issurjective() );
      assert( leer.isstrictlyascending() );
      
      assert(   leer.comparablewith( leer ) );
      assert(   leer.equals( leer ) );
      assert( ! leer.less( leer ) );
      
    }
    
    {

        cout << "Test Injection and Surjection" << endl;
        
        cout << "Injection" << endl;
        
        IndexMap inj( irB, irA );
        inj[3] = 2; inj[4] = 3; inj[5] = 4;
        
        inj.check();
        assert( inj.isinjective() );
        assert( !inj.issurjective() );
        assert( inj.getSourceRange() == irB );
        assert( inj.getDestRange() == irA );
        assert( inj.rangecontains( 4 ) );
        assert( !inj.rangecontains( 5 ) );
        
        cout << "Surjection" << endl;
        
        IndexMap sur( irD, irB );
        sur[0] = 4; sur[1] = 3; sur[2] = 5; sur[3] = 4;
        sur[4] = 3; sur[5] = 5; sur[6] = 4; sur[7] = 3;

        sur.check();
        assert( !sur.isinjective() );
        assert( sur.issurjective() );
        assert( sur.getSourceRange() == irD );
        assert( sur.getDestRange() == irB );
        
        cout << "Test composition of injection with surjection" << endl;
        const IndexMap prod = inj * sur;
        cout << inj << sur << prod << endl;
        
        IndexMap test( irD, irA );
        test[0] = 3; test[1] = 2; test[2] = 4; test[3] = 3;
        test[4] = 2; test[5] = 4; test[6] = 3; test[7] = 2;
        
        assert( prod.comparablewith( test ) );
        assert( prod == test );
        
        assert( test.rangecontains( 3 ) );
        assert( test.rangeposition( 3 ) == 0 );
        assert( test.rangecontains( 2 ) );
        assert( test.rangeposition( 2 ) == 1 );
        assert( ! test.rangecontains( 5 ) );
        
    }

    cout << "Finished Unit Test" << endl;

    return 0;
}
