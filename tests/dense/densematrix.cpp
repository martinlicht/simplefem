

/**/

#include <ostream>
#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Dense Matrix class" << nl;

    DenseMatrix A( 3, 4 );
    A.randommatrix();

    LOG << A << nl;
    LOG << 3 * A << nl;

    DenseMatrix B( 3, 4 );
    B.randommatrix();

    LOG << B << nl;
    LOG << A + B << nl;

    DenseMatrix I3(3,3);
    I3.unitmatrix();
    DenseMatrix I4(4,4);
    I4.unitmatrix();
    LOG << I3 << I4 << nl;
    LOG << I3 * A << nl;
    LOG << A * I4 << nl;
    auto S5 = 5. * I3;
    LOG << S5 << nl;
    LOG << S5 * A << nl;

    LOG << "Finished Unit Test" << nl;

    return 0;
}
