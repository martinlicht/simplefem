

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/global.lagrangemassmatrix.hpp"
#include "../../fem/global.lagrangestiffnessmatrix.hpp"
#include "../../fem/global.lagrangebrokenmassmatrix.hpp"
#include "../../fem/global.lagrangebrokenstiffnessmatrix.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Lagrange Mass Matrix" << endl;
        
        MeshSimplicial2D M = UnitSquare2D();
        
        M.check();
        
        cout << "Refinement..." << endl;
        
        int number_of_refinements = 0;
        for( int i = 0; i < number_of_refinements; i++ )
            M.uniformrefinement();
        
        cout << "...done" << endl;
        
        M.check();
        
        SparseMatrix massmatrix      = LagrangeMassMatrix( M, 1 );
        SparseMatrix stiffnessmatrix = LagrangeStiffnessMatrix( M, 1 );
        
        SparseMatrix broken_massmatrix      = LagrangeBrokenMassMatrix( M, 1 );
        SparseMatrix broken_stiffnessmatrix = LagrangeBrokenStiffnessMatrix( M, 1 );
        
        cout << massmatrix << endl;
        cout << stiffnessmatrix << endl;
        
        cout << broken_massmatrix << endl;
        cout << broken_stiffnessmatrix << endl;
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
