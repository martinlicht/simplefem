

/**/

#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"
#include "generateindexmaps.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Index Map Generators" << endl;
	
	if(false){
		
		cout << "Teste generator for all" << endl;
		
		IndexRange bereich( 1,3 );
		std::vector<IndexMap>  all = generateIndexMaps( bereich, bereich );

		for( const IndexMap& im : all )
			cout << im << endl;
		
	}
	
	if(false){
		
		cout << "Teste generator for permutations" << endl;
		
		IndexRange bereich( 1,3 );
		std::vector<IndexMap>  all = generatePermutations( bereich );

		for( const IndexMap& im : all )
			cout << im << endl;
		
	}
	cout << "Finished Unit Test" << endl;

	{
		
		cout << "Teste generator for sigmas" << endl;
		
		IndexRange from( 1,3 );
		IndexRange targ( 2,5 );
		std::vector<IndexMap>  all = generateSigmas( from, targ );

		for( const IndexMap& im : all )
			cout << im << endl;
		
	}
	


	return 0;
}
