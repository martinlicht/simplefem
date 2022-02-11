

/**/

#include <iostream>
// #include <fstream>
// #include <iomanip>

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
        
        // LOG << std::setprecision(10);
        
        const int n_min = 1;
        const int n_max = 3;
        const int r_min = 1;
        const int r_max = 6;

        for( int n = n_min; n <= n_max; n++ )
        for( int r = r_min; r <= r_max; r++ )
        {
            LOG << "Dimension: " << space << n_min << " <= " << n << " <= " << n_max << endl;
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << endl;
            
            const auto lpsbc = InterpolationPointsBarycentricCoordinates( n, r );
            
            const auto EM = EvaluationMatrix( r, lpsbc );
            
            assert( EM.issquare() );
        
            const auto EMinv = Inverse( EM );
        
            int N = EM.getdimin();

            Float diff_inv_1 = ( EM * EMinv - IdentityMatrix(N) ).norm();
            Float diff_inv_2 = ( EMinv * EM - IdentityMatrix(N) ).norm();
            
            LOG << "dim(A)=" << N << "\tdiff1=" << diff_inv_1 << "\tdiff2=" << diff_inv_2 << space << machine_epsilon << nl;
            
            assert( diff_inv_1 < std::sqrt( machine_epsilon ) );
            assert( diff_inv_2 < std::sqrt( machine_epsilon ) );
            
        }
        
        LOG << "Finished Unit Test" << endl;
        
        
        return 0;
}
