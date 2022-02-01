

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../dense/cholesky.hpp"
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

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Inverse of Poly Matrix" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(10);

        const int r_min = 1;
        const int r_max = 9;

        const int n_min = 1;
        const int n_max = 3;

        for( int n = n_min; n <= n_max; n++ )
        for( int r = r_min; r <  r_max; r++ )
        {
            
            DenseMatrix MM = polynomialmassmatrix( n, r );
        
            int N = MM.getdimin();

            LOG << "Dimension: " << space << n_min << " <= " << n << " <= " << n_max << endl;
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << endl;

                    
            LOG << "Matrix dimension: " << N << nl;
            LOG << "Determinant: " << Determinant(MM) << nl;

            LOG << "Inverse..." << nl;
            DenseMatrix MMinv = Inverse(MM);

            LOG << "Cholesky decomposition..." << nl;
            DenseMatrix MMchol = CholeskyDecomposition(MM);
            
            LOG << "QR decomposition..." << nl;
            DenseMatrix MMqr_q(MM), MMqr_r(MM);
            QRFactorization(MM,MMqr_q,MMqr_r);
            
            Float diff_inv  = ( MM * MMinv - IdentityMatrix(N) ).norm();
            Float diff_chol = ( MMchol * Transpose(MMchol) - MM ).norm();
            Float diff_qr   = ( MMqr_q * MMqr_r - MM ).norm();

            LOG << "\ta=" << diff_inv
                << "\tb=" << diff_chol
                << "\tc=" << diff_qr
                << nl;
                        
        }
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
