

/**/

#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Index Ranges" << endl;
	
	cout << "Test invalid index ranges" << std::endl;
	
// 	assert( false );
	
	if( false ) {
	  IndexRange dummy( 0, std::numeric_limits<int>::max() );
	  cout << dummy << " " << std::numeric_limits<int>::max();
	  cout << "Invalid range created" << std::endl;
	}
	  
	cout << "Test empty index ranges" << std::endl;
	
	IndexRange irE( 3, 2 );
	
	assert( irE.isempty() );
	assert( irE.cardinality() == 0 );
	assert( irE.cardinality() == irE.getlength() );
	
	cout << "Test non-empty index ranges" << std::endl;
	
	IndexRange irA( 3, 7 );
	IndexRange irB( 5, 5 );
	
	assert( irA.cardinality() == 5 );
	
	assert( irA.getlength() == irA.cardinality() );
	
	assert( !irA.contains(2) );
	assert( irA.contains(3) );
	assert( irA.contains(4) );
	assert( irA.contains(5) );
	assert( irA.contains(6) );
	assert( irA.contains(7) );
	assert( !irA.contains(8) );
	
	assert( irB.getlength() == irB.cardinality() );
	
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
	
// 	cout << "Test combination of index ranges" << std::endl;
// 	
// 	
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
