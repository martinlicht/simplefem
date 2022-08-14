#include <ostream>
// #include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.whitneyincl.hpp"

#include "../../fem/lagrangematrices.hpp"



using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "assemble some common matrices" );

int main()
{
    LOG << "Unit Test: " << TestName << endl;
        
    {
        LOG << "Case 1D" << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        LOG << "Refinement..." << endl;
        
        int number_of_refinements = 2;
        
        for( int i = 0; i < number_of_refinements; i++ )
            M.uniformrefinement();
        
        LOG << "...assemble matrices" << endl;
        
        for( int r = 0; r <= 4; r++ ) 
        for( int k = 0; k <= 1; k++ ) 
        {
            LOG << "[ k, r ] = [" << k << ", " << r << "]\n";
            
            SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
            
            SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), k, r );
            
            if( r >= 1 && k < M.getinnerdimension() )
                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r );
            
            if( r == 1 && k == 0 )
                SparseMatrix inclmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );
            
            for( int r_plus = 0; r_plus < 5; r_plus++ )
                SparseMatrix elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r, r_plus );
            
            if( r > 0 )
                SparseMatrix SullivanInclMatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            if( r > 0 )
                SparseMatrix WhitneyInclMatrix = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            SparseMatrix lagrange_massmatrix      = LagrangeMassMatrix( M, 1 );
            SparseMatrix lagrange_stiffnessmatrix = LagrangeStiffnessMatrix( M, 1 );
            
            SparseMatrix lagrange_broken_massmatrix      = LagrangeBrokenMassMatrix( M, 1 );
            SparseMatrix lagrange_broken_stiffnessmatrix = LagrangeBrokenStiffnessMatrix( M, 1 );
            
        }
    }

    {
        LOG << "Case 2D" << endl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        LOG << "Refinement..." << endl;
        
        int number_of_refinements = 3;
        
        for( int i = 0; i < number_of_refinements; i++ )
            M.uniformrefinement();
        
        LOG << "...assemble matrices" << endl;
        
        for( int r = 0; r <= 3; r++ ) 
        for( int k = 0; k <= 2; k++ ) 
        {
            LOG << "[ k, r ] = [" << k << ", " << r << "]\n";
            
            SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
            
            SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), k, r );
            
            if( r >= 1 && k < M.getinnerdimension() )
                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r );
            
            if( r == 1 && k == 0 )
                SparseMatrix inclmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );
            
            for( int r_plus = 0; r_plus < 5; r_plus++ )
                SparseMatrix elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r, r_plus );
            
            if( r > 0 )
                SparseMatrix SullivanInclMatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            if( r > 0 )
                SparseMatrix WhitneyInclMatrix = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), k, r );

            SparseMatrix lagrange_massmatrix      = LagrangeMassMatrix( M, 1 );
            SparseMatrix lagrange_stiffnessmatrix = LagrangeStiffnessMatrix( M, 1 );
            
            SparseMatrix lagrange_broken_massmatrix      = LagrangeBrokenMassMatrix( M, 1 );
            SparseMatrix lagrange_broken_stiffnessmatrix = LagrangeBrokenStiffnessMatrix( M, 1 );
            
        }
    }

    {
        LOG << "Case 3D" << endl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        LOG << "Refinement..." << endl;
        
        int number_of_refinements = 1;
        
        for( int i = 0; i <= number_of_refinements; i++ ) 
            M.uniformrefinement();
        
        LOG << "...assemble matrices" << endl;
        
        for( int r = 0; r <= 3; r++ ) 
        for( int k = 0; k <= 3; k++ ) 
        {
            LOG << "[ k, r ] = [" << k << ", " << r << "]\n";
            
            SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
            
            SparseMatrix massmatrix_rhs = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), k, r );
            
            if( r >= 1 && k < M.getinnerdimension() )
                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r );
            
            if( r == 1 && k == 0 )
                SparseMatrix inclmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );
            
            for( int r_plus = 0; r_plus < 5; r_plus++ )
                SparseMatrix elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r, r_plus );
            
            if( r > 0 )
                SparseMatrix SullivanInclMatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );

            if( r > 0 )
                SparseMatrix WhitneyInclMatrix = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), k, r );

            SparseMatrix lagrange_massmatrix      = LagrangeMassMatrix( M, 1 );
            SparseMatrix lagrange_stiffnessmatrix = LagrangeStiffnessMatrix( M, 1 );
            
            SparseMatrix lagrange_broken_massmatrix      = LagrangeBrokenMassMatrix( M, 1 );
            SparseMatrix lagrange_broken_stiffnessmatrix = LagrangeBrokenStiffnessMatrix( M, 1 );
            
        }
    }
        
    LOG << "Finished Unit Test: " << TestName << endl;
    
    return 0;
}
