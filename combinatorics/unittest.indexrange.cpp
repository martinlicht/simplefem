

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
	
	try {
	  IndexRange dummy( 0, std::numeric_limits<int>::max() );
	} catch (...) {
	  (void)0;
	}
	
	cout << "Test empty index ranges" << std::endl;
	
	IndexRange irE( 3, 2 );
	
	attest( irE.isempty() );
	attest( irE.cardinality() == 0 );
	attest( irE.cardinality() == irE.getlength() );
	
	cout << "Test non-empty index ranges" << std::endl;
	
	IndexRange irA( 3, 7 );
	IndexRange irB( 5, 5 );
	
	attest( irA.cardinality() == 5 );
	
	attest( irA.getlength() == irA.cardinality() );
	
	attest( !irA.contains(2) );
	attest( irA.contains(3) );
	attest( irA.contains(4) );
	attest( irA.contains(5) );
	attest( irA.contains(6) );
	attest( irA.contains(7) );
	attest( !irA.contains(8) );
	
	attest( irB.getlength() == irB.cardinality() );
	
	attest( irB.cardinality() == 1 );
	
	attest( irB.contains(5) );
	attest( !irB.contains(4) );
	attest( !irB.contains(6) );
	
	cout << "Test indexing in non-empty index ranges" << std::endl;
	
	attest( irB.element2position(5) == 0 );
	attest( irB.position2element(0) == 5 );
	attest( irA.element2position(3) == 0 );
	attest( irA.element2position(5) == 2 );
	attest( irA.element2position(7) == 4 );
	attest( irA.position2element(0) == 3 );
	attest( irA.position2element(2) == 5 );
	attest( irA.position2element(4) == 7 );
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
