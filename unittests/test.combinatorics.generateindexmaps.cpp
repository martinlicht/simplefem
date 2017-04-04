

/**/

#include <iostream>
#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Index Map Generators" << endl;
        
        if(true){
                
                cout << "Teste generator for empty index map" << endl;
                
                IndexRange from( 1,-3 );
                IndexRange targ( 2,5 );
                
                std::vector<IndexMap>  all = generateEmptyMap( from, targ );
                
                for( const IndexMap& im : all )
                        cout << im << endl;
                
                cout << "Tested" << endl;
                
        }
        
        if(true){
                
                cout << "Teste generator for all: [1,3] -> [1,3]" << endl;
                
                IndexRange bereich( 1,3 );
                std::vector<IndexMap> all = generateIndexMaps( bereich, bereich );

                for( const IndexMap& im : all )
                        cout << im << endl;
                
                cout << "Tested" << endl;
                
        }
        
        if(true){
                
                cout << "Teste generator for permutations of [1,3]" << endl;
                
                IndexRange bereich( 1,3 );
                std::vector<IndexMap>  all = generatePermutations( bereich );

                for( const IndexMap& im : all )
                        cout << im << endl << "signum: " << signPermutation( im ) << endl;
                
                cout << "Tested" << endl;
                
        }
        
        if(true){
                
                cout << "Teste generator for sigmas [1,3] -> [2,5]" << endl;
                
                IndexRange from( 1,3 );
                IndexRange targ( 2,5 );
                std::vector<IndexMap>  all = generateSigmas( from, targ );

                for( const IndexMap& im : all )
                        cout << im << endl;
                
                cout << "Tested" << endl;
                
        }
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
