

/**/

#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"
#include "multiindex.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Multi-Indices" << endl;
	
	IndexRange irA( 2, 5 );
	MultiIndex miA( irA );
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
