

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
        cout << "Unit Test for Inverse of Evaluation Matrix" << endl;
        
        cout << std::setprecision(10);

        for( int n = 1; n <= 3; n++ )
        for( int r = 1; r <  7; r++ )
        {
            
            const auto lpsbc = InterpolationPointsBarycentricCoordinates( n, r );
            
            const auto EM = EvaluationMatrix( n, r, lpsbc );
        
            const auto EMinv = Inverse( EM );
        
            int N = EM.getdimin();

            Float diff_inv = ( EM * EMinv - IdentityMatrix(N) ).norm();//
            
            std::cout << "\tn=" << n
                      << "\tr=" << r
                      << "\tN=" << N
                      << "\ta=" << diff_inv
                      << nl;

        }
        
        cout << "Finished Unit Test" << endl;
        
        
        return 0;
}
