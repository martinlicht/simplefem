

/**/

#include <iostream>
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"
#include "../../combinatorics/multiindex.hpp"


using namespace std;

int main()
{
    
    LOG << "Unit Test for Multi-Indices";// << endl;
        
    if( true )
    {
        LOG << "Basic functionality";// << endl;
        
        IndexRange irA( 2, 5 );
        MultiIndex miA( irA );
        
        assert( miA.absolute() == 0 );
        assert( miA.factorial() == 1 );
        
        miA += 4;

        assert( miA.absolute() == 1 );
        assert( miA.factorial() == 1 );
        
        miA += 4;
        
        assert( miA.absolute() == 2 );
        assert( miA.factorial() == 2 );
        
        miA += 2;

        assert( miA.absolute() == 3 );
        assert( miA.factorial() == 2 );
        
        miA += 2;
        
        assert( miA.absolute() == 4 );
        assert( miA.factorial() == 4 );
        
        miA += 2;

        assert( miA[2] == 3 and miA[3] == 0 and miA[4] == 2 and miA[5] == 0 );
        assert( miA.absolute() == 5 );
        assert( miA.factorial() == 2 * 6 );
        
        miA[3] = 5;
        
        assert( miA[2] == 3 and miA[3] == 5 and miA[4] == 2 and miA[5] == 0 );
        assert( miA.absolute() == 10 );
        assert( miA.factorial() == 2 * 6 * 5*4*3*2*1 );
        
    }
    
    if( true )
    {
        LOG << "Comparison and arithmetics";// << endl;
        
        IndexRange irA( 2, 5 );
        IndexRange irB( 1, 4 );
        
        MultiIndex miA1( irA );
        MultiIndex miA2( irA );
        MultiIndex miB ( irB );
        
        assert( miA1.comparablewith( miA2 ) );
        assert( miA2.comparablewith( miA1 ) );
        assert( not miA1.comparablewith( miB ) );
        assert( not miB.comparablewith( miA1 ) );
        
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
        
        LOG << "Comparison and arithmetics done";// << endl;
    }
    
    LOG << "Finished Unit Test";// << endl;
    
    return 0;
}
