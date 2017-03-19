

/**/

#include <iostream>
#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Multi-Indices" << endl;
	
	IndexRange irA( 2, 5 );
	MultiIndex miA( irA );
	
	cout << miA << endl;
	miA += 4;
	miA += 2;
	miA += 2;
	miA += 2;
	cout << miA << endl;
	cout << miA.absolute() << space << miA.factorial() << endl;
	
	miA[3] = 7;
	cout << miA << endl;
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
