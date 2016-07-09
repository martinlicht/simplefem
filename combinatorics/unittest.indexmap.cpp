

/**/

#include <iostream>
#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for Index Mapping" << endl;
	
	{
		
		cout << "Teste IndexRange class" << endl;
		
		IndexRange einmal( 2, 5 );
		IndexRange zweimal( 3, 5 );
		
		cout << einmal << zweimal << endl;
		cout << einmal.contains(2) << " " << zweimal.contains(2) << endl;
		cout << einmal.contains( zweimal ) << endl;
		
	}
	
	{
	
		cout << "Teste IndexMappings" << endl;
		
		IndexRange vorne( 2,3 );
		IndexRange hinten( 0,4 );
		
		cout << "Abbildungen" << endl;
		IndexMap inj( vorne, hinten );
		IndexMap sur( hinten, vorne );
		cout << "daten laden..." << endl;
		inj[2] = 0; inj[3] = 4;
		sur[0] = 2; sur[1] = 3; sur[2] = 2; sur[3] = 2; sur[4] = 3;
		
		cout << "checken..." << endl;
		inj.check();
		sur.check();
		
		cout << "eigenschaften..." << endl;
		cout << inj << endl << sur << endl;
		cout << inj.isinjective() << " " << inj.issurjective() << endl;
		cout << sur.isinjective() << " " << sur.issurjective() << endl;
		
		cout << "product..." << endl;
		IndexMap prod = sur * inj;
		cout << prod << endl;
		cout << prod.isinjective() << " " << prod.issurjective() << endl;
		cout << prod.skip(3) << endl;
		
	}
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
