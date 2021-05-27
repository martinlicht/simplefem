

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test: Evaluation Matrix and its Invertibility" << endl;
        
        LOG << std::setprecision(10);

        for( int n = 1; n <= 3; n++ )
        for( int r = 1; r <  7; r++ )
        {
            
            const auto lpsbc = InterpolationPointsBarycentricCoordinates( n, r );
            
            const auto EM = EvaluationMatrix( r, lpsbc );
            
            assert( EM.issquare() );
        
            const auto EMinv = Inverse( EM );
        
            int N = EM.getdimin();

            Float diff_inv_1 = ( EM * EMinv - IdentityMatrix(N) ).norm();
            Float diff_inv_2 = ( EMinv * EM - IdentityMatrix(N) ).norm();
            
            LOG << "\tn=" << n
                      << "\tr=" << r
                      << "\tN=" << N
                      << "\tdiff1=" << diff_inv_1
                      << "\tdiff2=" << diff_inv_2
                      << nl;

        }
        
        LOG << "Finished Unit Test" << endl;
        
        
        return 0;
}
