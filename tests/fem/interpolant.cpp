

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        LOG << "Unit Test for Interpolation in FEEC" << endl;
        
        MeshSimplicial2D M = UnitDisk(6);
        
        M.check();
        
        auto scalarfield = [](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            return FloatVector({ std::sqrt( vec[0]*vec[0] + vec[1]*vec[1] ) });
            return FloatVector({ 1 + vec[0] + vec[1] * vec[1] });
        };
        
        FloatVector results = Interpolation( M, M.getinnerdimension(), 0, 2, scalarfield );
        
        LOG << results << endl;
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
