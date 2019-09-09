

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"
#include "../../combinatorics/multiindex.hpp"


using namespace std;

int main()
{
    
    cout << "Unit Test for Multi-Indices" << endl;
        
    {
        cout << "Basic functionality" << endl;
        
        IndexRange irA( 2, 5 );
        MultiIndex miA( irA );
        
        cout << "Before adding" << endl;
        cout << miA << endl;
        miA += 4;
        miA += 4;
        miA += 2;
        miA += 2;
        miA += 2;
        cout << "After adding" << endl;
        cout << miA << endl;
        cout << "Absolute and factorial" << endl;
        cout << miA.absolute() << space << miA.factorial() << endl;
        cout << endl;
        
        cout << "Assignment" << endl;
        miA[3] = 7;
        cout << miA << endl;
    }
    
    {
        cout << "Comparison and arithmetics" << endl;
        
        IndexRange irA( 2, 5 );
        IndexRange irB( 1, 4 );
        MultiIndex miA1( irA );
        MultiIndex miA2( irA );
        MultiIndex miB ( irB );
        
        assert( miA1.comparablewith( miA2 ) );
        assert( miA2.comparablewith( miA1 ) );
        assert( ! miA1.comparablewith( miB ) );
        assert( ! miB.comparablewith( miA1 ) );
        
        miA1 += 4;
        miA1 += 4;
        miA1 += 2;
        miA1 += 2;
        miA1 += 2;
        
        miA2 += 2;
        miA2 += 4;
        miA2 += 2;
        miA2 += 4;
        miA2 += 2;
        
        assert( miA1.comparablewith( miA2 ) );
        assert( miA2.comparablewith( miA1 ) );
        assert( miA1 == miA2 );
        assert( miA2 == miA1 );
        
        miA2 -= 4;
        
        assert( miA1.comparablewith( miA2 ) );
        assert( miA2.comparablewith( miA1 ) );
        assert( miA1 != miA2 );
        assert( miA2 != miA1 );
        
        miA1 -= 4;
        
        assert( miA1.comparablewith( miA2 ) );
        assert( miA2.comparablewith( miA1 ) );
        assert( miA1 == miA2 );
        assert( miA2 == miA1 );
        
        MultiIndex miA3( irA );
        miA3[2] = 6;
        miA3[3] = 0;
        miA3[4] = 2;
        miA3[5] = 0;

        assert( miA3.comparablewith( miA1+miA2 ) );
        assert( ( miA1+miA2 ).comparablewith( miA3 ) );
        assert( miA3 == ( miA1+miA2 ) );
        assert( ( miA1+miA2 ) == miA3 );
        
        assert( ( miA1-miA2 ).at(2) == 0 );
        assert( ( miA1-miA2 ).at(3) == 0 );
        assert( ( miA1-miA2 ).at(4) == 0 );
        assert( ( miA1-miA2 ).at(5) == 0 );
        
        cout << "Comparison and arithmetics done" << endl;
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
