

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../dense/qr.factorization.hpp"
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
        cout << "Unit Test for Inverse of Poly Matrix" << endl;
        
        cout << std::setprecision(10);

        for( int n = 1; n <= 3; n++ )
        for( int r = 1; r <  10; r++ )
        {
            
            DenseMatrix MM = polynomialmassmatrix( n, r );
        
            DenseMatrix MMinv = Inverse(MM);

            DenseMatrix MMchol = CholeskyDecomposition(MM);
            
            DenseMatrix MMqr_q(MM), MMqr_r(MM);
            QRFactorization(MM,MMqr_q,MMqr_r);
            
            int N = MM.getdimin();

            Float diff_inv  = ( MM * MMinv - IdentityMatrix(N) ).norm();//
            Float diff_chol = ( MMchol * Transpose(MMchol) - MM ).norm();
            Float diff_qr   = ( MMqr_q * MMqr_r - MM ).norm();

            std::cout << "\tn=" << n
                      << "\tr=" << r
                      << "\tN=" << N
                      << "\ta=" << diff_inv
                      << "\tb=" << diff_chol
                      << "\tc=" << diff_qr
                      << nl;
                      
            // if( n==2 and r==2 )
                // std::cout << MM / factorial_integer(2) << nl;;

        }
        
        cout << "Finished Unit Test" << endl;
        
        
        
        return 0;
}
