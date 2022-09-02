

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../sparse/sparsematrix.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for SparseMatrix" << nl;

    SparseMatrix M( 2, 3 );

    for( int i = 0; i < 5; i++ )
        for( int j = 0; j < 7; j++ )
            M.addentry( (3*i) % 2, (2*j) % 3, i / 3. + j*j );

    LOG << "This is the content of some matrix:" << nl;
    LOG << M << nl;

    LOG << "Sort the Entries" << nl;
    M.sortentries();
    LOG << M << nl;

    M.clearentries();
    LOG << "Empty Matrix again" << nl;
    LOG << M << nl;

    LOG << "Next Matrix:" << nl;
    M.addentry( 0, 0, 1. );
    M.addentry( 0, 1, 2. );
    M.addentry( 0, 2, 3. );
    M.addentry( 1, 0, 4. );
    M.addentry( 1, 1, 5. );
    M.addentry( 1, 2, 6. );
    LOG << M << nl;

    FloatVector vec(3);
    vec.setentry(0,13);
    vec.setentry(1,17);
    vec.setentry(2,19);
    LOG << "Some vector:" << nl;
    LOG << vec << nl;

    LOG << "Matrix-Vector Product:" << nl;
    LOG << M * vec << nl;

    LOG << "Finished Unit Test" << nl;

    return 0;
}
