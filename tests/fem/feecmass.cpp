

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/global.sullivanincl.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for FEEC Mass Matrix" << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        cout << "Refinement..." << endl;
        
        int number_of_refinements = 3;
        for( int i = 0; i < number_of_refinements; i++ )
            M.uniformrefinement();
        
        cout << "...done" << endl;
        
        cout << "...assemble matrices" << endl;
        
        for( int r = 0; r <= 5; r++ ) 
        for( int k = 0; k <= 2; k++ ) 
        {
            cout << "[ k, r ] = [" << k << ", " << r << "]\n";
            
            SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
            
            SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), k, r );
            
            if( r >= 1 && k < M.getinnerdimension() )
            SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r );
            
            if( r == 1 && k == 0 )
            SparseMatrix inclmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );
            
            SparseMatrix elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, 4, 3 );
            
            if( r > 0 )
                SparseMatrix SullivanInclMatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
        }
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
