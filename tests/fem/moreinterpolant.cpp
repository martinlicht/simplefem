

/**/

#include <iostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Interpolation in FEEC" << endl;
        
        cout << "2D Calculations" << endl;
        
        {

            MeshSimplicial2D M = UnitDisk(3);
            
            cout << "... mesh done" << endl;
        
            auto scalarfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ sqrt( vec[0]*vec[0] + vec[1]*vec[1] ) });
                return FloatVector({ 1 + vec[0] + vec[1] * vec[1] });
            };
            
            auto vectorfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 1 + vec[0], vec[1] * vec[1] });
            };
            
            cout << "\n... k=0" << endl;
            for( int r = 0; r <  3; r++ )
                Interpolation( M, M.getinnerdimension(), 0, r, scalarfield ),
                cout << " r=" << r;
            
            cout << "\n... k=1" << endl;
            for( int r = 0; r <  3; r++ )
                Interpolation( M, M.getinnerdimension(), 1, r, vectorfield ),
                cout << " r=" << r;
            
            cout << "\n... k=2" << endl;
            for( int r = 0; r <  3; r++ )
                Interpolation( M, M.getinnerdimension(), 2, r, scalarfield ),
                cout << " r=" << r;

        }
        
        cout << "\n3D Calculations" << endl;
        
        {

            MeshSimplicial3D M = UnitSimplex3D();
//             M.uniformrefinement();
            
            cout << "... mesh done" << endl;
            
            auto scalarfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ sqrt( vec[0]*vec[0] + vec[1]*vec[2] ) });
            };
            
            auto vectorfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ 1 + vec[0], vec[1] * vec[1], 2. * vec[2] });
            };
            
            int Rmax = 2;
            
            cout << "... k=0" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );
            
            cout << "... k=1" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );
            
            cout << "... k=2" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 2, r, vectorfield );

            cout << "... k=3" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 3, r, scalarfield );

        }
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
