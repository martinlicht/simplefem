

/**/

#include <iostream>
#include "../basic.hpp"
#include "floatvector.hpp"
#include "sparsematrix.hpp"


using namespace std;

int main()
{
	cout << "Unit Test for SparseMatrix" << endl;
	
	SparseMatrix M( 2, 3 );
	
	for( int i = 0; i < 5; i++ )
		for( int j = 0; j < 7; j++ )
			M.addentry( (3*i) % 2, (2*j) % 3, i / 3. + j*j );
	
	cout << "This is the content of some matrix:" << endl;
	cout << M << endl;
	
	cout << "Sort the Entries" << endl;
	M.sortentries();
	cout << M << endl;
	
	M.clearentries();
	cout << "Empty Matrix again" << endl;
	cout << M << endl;
	
	cout << "Next Matrix:" << endl;
	M.addentry( 0, 0, 1. );
	M.addentry( 0, 1, 2. );
	M.addentry( 0, 2, 3. );
	M.addentry( 1, 0, 4. );
	M.addentry( 1, 1, 5. );
	M.addentry( 1, 2, 6. );
	cout << M << endl;
	
	FloatVector vec(3);
	vec.setentry(0,13);
	vec.setentry(1,17);
	vec.setentry(2,19);
	cout << "Some vector:" << endl;
	cout << vec << endl;
	
	cout << "Matrix-Vector Product:" << endl;
	cout << M * vec << endl;
	
	cout << "Finished Unit Test" << endl;

	return 0;
}
