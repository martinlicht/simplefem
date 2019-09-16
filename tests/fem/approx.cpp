

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
#include "../../fem/foo.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Inverse of Poly Matrix" << endl;
        
        cout << std::setprecision(10);

        for( int n = 1; n <= 3; n++ )
        for( int r = 1; r <  7; r++ )
        {
            
            DenseMatrix MM = polynomialmassmatrix( n, r );
        
            DenseMatrix MMinv = Inverse(MM);

            DenseMatrix MMfac = CholeskyDecomposition(MM);
            
            int N = MM.getdimin();

            Float diff_inv = ( MM * MMinv - IdentityMatrix(N) ).norm();//
            Float diff_fac = ( MMfac * Transpose(MMfac) - IdentityMatrix(N) ).norm();

            auto allone = FloatVector( N, 1. );

            auto test = allone.scalarproductwith( MM * allone );

            std::cout << "\tn=" << n
                      << "\tr=" << r
                      << "\tN=" << N
                      << "\ta=" << diff_inv
                      << "\tb=" << diff_fac
                      << "\tc=" << 1./test << ' ' << test
                      << nl;
                      
            if( n==2 and r==2 )
                std::cout << MM / factorial_integer(2) << nl;;

        }
        
        cout << "Finished Unit Test" << endl;
        
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
        
        //return 0;
        
        
        
        cout << "Unit Test for Approximating L2 norm of scalar field" << endl;
        
        cout << std::setprecision(10);

        if(false){

            cout << "Case 2D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = UnitSquare2D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_field;
            std::vector<Float>                                          experiments_value;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
                }
            );

            experiments_value.push_back(3.194528049465325113615213730287503906590157785275923662043);
            
            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ ( vec[0] * vec[1] > 0 ) ? std::exp( vec[0] ) : 0. });
                }
            );

            experiments_value.push_back(3.626860407847018767668213982801261704886342012321135721309);
            

            
            cout << "Calculate L2 products" << endl;
            
            for( int l = 0; l <= 3; l++ ){
                
                cout << "Refinement..." << endl;
            
                M.uniformrefinement();
                
                cout << "...assemble matrices" << endl;
            
                for( int r = 0; r <= 8; r++ ) 
                {
                    SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                    
                    for( int i = 0; i < experiments_field.size(); i++){

                        const auto& scalarfield = experiments_field[i];
                        const auto& should_be   = experiments_value[i];

                        FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                        Float mass1 = interpol.scalarproductwith( massmatrix * interpol );
                        //Float mass2 = power( ( massmatrix_rhs * interpol ).norm(), 2. );

                        cout << "[ i, l, r ] = [" << i << ", "  << l << ", " << r << "]\t";
                        
                        cout << "mass: " << mass1 << ' ' << ' ' << should_be << endl;// << mass2

                    }
                    
                }
            } 
        
        }
        
        
        {

            cout << "Case 3D" << endl;
            
            cout << "Initial mesh..." << endl;
            
            MeshSimplicial3D M = UnitSimplex3D();
            
            M.check();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_field;
            std::vector<Float>                                          experiments_value;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            experiments_field.push_back( 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({ 1. });
                }
            );

            experiments_value.push_back(1.);
            
            
            

            
            cout << "Calculate L2 products" << endl;
            
            for( int l = 0; l <= 2; l++ ){
                
                cout << "Refinement..." << endl;
            
                M.uniformrefinement();
                
                cout << "...assemble matrices" << endl;
            
                for( int r = 0; r <= 7; r++ ) 
                {
                    SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    //SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                    
                    for( int i = 0; i < experiments_field.size(); i++){

                        const auto& scalarfield = experiments_field[i];
                        const auto& should_be   = experiments_value[i];

                        FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                        Float mass1 = interpol.scalarproductwith( massmatrix * interpol );
                        Float mass2 = 0.;//power( ( massmatrix_rhs * interpol ).norm(), 2. );

                        cout << "[ i, l, r ] = [" << i << ", "  << l << ", " << r << "]\t";
                        
                        cout << "mass: " << mass1 << ' ' << mass2 << ' ' << should_be << endl;

                    }
                    
                }
            } 
        
        }
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
