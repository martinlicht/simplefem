

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/densematrix.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: (3D) degree elevation commutes with exterior derivative";// << endl;
        
        LOG << std::setprecision(10);

        LOG << "Initial mesh...";// << endl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> fields;
        
        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( vec[0] + vec[1] - vec[2] ) });
            }
        );

        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                auto ret = FloatVector({ 
                                std::exp( vec[0] + vec[2] ), 
                                std::sin( -5.*vec[0] -vec[1] ),
                                std::atan( vec[1] * vec[2] )
                });
                assert( ret.getdimension() == 3 );
                return ret;
            }
        );

        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                auto ret = FloatVector({ 
                                std::cos( -5.*vec[0]*vec[0] ), 
                               -std::exp( 5.*vec[2] / ( 1. + vec[0] * vec[0] ) ), 
                                std::cos( 2.*vec[1]-vec[0] ) 
                                
                });
                assert( ret.getdimension() == 3 );
                return ret;
            }
        );
        
        

        const int r_min = 1;
        
        const int r_max = 3;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        const int r_plus_max = 3;
         
        Float errors[ M.getinnerdimension() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            for( int k      =     0; k      <  M.getinnerdimension(); k++      ) 
            for( int r      = r_min; r      <=                 r_max; r++      ) 
            for( int r_plus =     0; r_plus <=            r_plus_max; r_plus++ ) 
            {
                
                LOG << "...assemble matrices: l=" << l << " k=" << k << " r=" << r << " rplus=" << r_plus;// << endl;
        
                SparseMatrix lower_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r          );

                SparseMatrix upper_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r + r_plus );

                SparseMatrix diyi_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k  , r  , r_plus );
                
                SparseMatrix dier_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k+1, r-1, r_plus );
                
                SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k+1, r + r_plus - 1 );
                
                
                
                FloatVector interpol_function = Interpolation( M, M.getinnerdimension(), k, r, fields[k] );

                auto path1 = dier_elevation * lower_diffmatrix * interpol_function;

                auto path2 = upper_diffmatrix * diyi_elevation * interpol_function;

                auto commutator_error = path1 - path2;
                
                Float commutator_error_mass = commutator_error * ( massmatrix * commutator_error );

                assert( std::isfinite( commutator_error_mass ) );
                
                errors[k][l-l_min][r-r_min][r_plus] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                
                
            }

            LOG << "Refinement...";// << endl;
        
            M.uniformrefinement();
            
            

        } 
        
        
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable[ M.getinnerdimension() ];
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < M.getinnerdimension(); i++ ) 
                    contable[i] << errors[i][l-l_min][r-r_min][r_plus_max];
                        
            }
            
            for( int i = 0; i < M.getinnerdimension(); i++ ) 
                contable[i] << nl; 
            
        }
        
        
        
        for( int i = 0; i < M.getinnerdimension(); i++ ) 
        {
            contable[i].lg(); 
            LOG << "-------------------" << nl;
        }
                
        
        
        
        
        LOG << "Check that differences are small" << nl;
        
        for( int l      = l_min; l      <=           l_max; l++      ) 
        for( int r      = r_min; r      <=           r_max; r++      ) 
        for( int r_plus =     0; r_plus <=      r_plus_max; r_plus++ ) 
        for( int i      =     0; i < M.getinnerdimension(); i++      ) 
        {
            assert( errors[i][l-l_min][r-r_min][r_plus] < 10e-6 );
        }
            
        
        LOG << "Finished Unit Test";// << endl;
        
        return 0;
}
