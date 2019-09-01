

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/local.massmatrix.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for FEEC Mass Matrix" << endl;
        
        MeshSimplicial2D M = UnitSquare2D();
        
        M.check();
        
        cout << "Refinement..." << endl;
        
        int number_of_refinements = 0;
        for( int i = 0; i < number_of_refinements; i++ )
            M.uniformrefinement();
        
        cout << "...done" << endl;
        
        M.check();
        
        cout << M.getcoordinates() << endl;
        
        SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, 1 );
        
//         cout << massmatrix << endl;
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
