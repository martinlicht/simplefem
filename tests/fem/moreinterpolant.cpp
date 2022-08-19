

/**/

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Interpolation in FEEC" << endl;
        
        LOG << "2D Calculations" << endl;
        
        {

            MeshSimplicial2D M = UnitDisk(3);
            
            LOG << "... mesh done" << endl;
        
            auto scalarfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ sqrt( vec[0]*vec[0] + vec[1]*vec[1] ) });
                return FloatVector({ 1 + vec[0] + vec[1] * vec[1] });
            };
            
            auto vectorfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 1 + vec[0], vec[1] * vec[1] });
            };
            
            LOG << "\n... k=0" << endl;
            for( int r = 0; r <  3; r++ ) {
                Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );
                LOG << " r=" << r;
            }
            
            LOG << "\n... k=1" << endl;
            for( int r = 0; r <  3; r++ ) {
                Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );
                LOG << " r=" << r;
            }
            
            LOG << "\n... k=2" << endl;
            for( int r = 0; r <  3; r++ ) {
                Interpolation( M, M.getinnerdimension(), 2, r, scalarfield );
                LOG << " r=" << r;
            }

        }
        
        LOG << "\n3D Calculations" << endl;
        
        {

            MeshSimplicial3D M = UnitSimplex3D();
//             M.uniformrefinement();
            
            LOG << "... mesh done" << endl;
            
            auto scalarfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ sqrt( vec[0]*vec[0] + vec[1]*vec[2] ) });
            };
            
            auto vectorfield = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ 1 + vec[0], vec[1] * vec[1], 2. * vec[2] });
            };
            
            int Rmax = 2;
            
            LOG << "... k=0" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );
            
            LOG << "... k=1" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );
            
            LOG << "... k=2" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 2, r, vectorfield );

            LOG << "... k=3" << endl;
        
            for( int r = 0; r <= Rmax; r++ )
                FloatVector results = Interpolation( M, M.getinnerdimension(), 3, r, scalarfield );

        }
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
