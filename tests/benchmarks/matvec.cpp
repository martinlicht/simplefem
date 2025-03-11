

/**/
#include <ctime>

#include "../basic.hpp"
#include "../mesh/coordinates.hpp"
#include "../mesh/mesh.simplicial2D.hpp"
#include "../mesh/examples2D.hpp"
#include "../sparse/matcsr.hpp"
#include "../fem/global.massmatrix.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Benchmark: (2D) CSR-matrix vector product" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D M = StandardSquare2D_simple();
    
    M.check();

    
    
    
    const int r = 4;
    
    const int l = 6;
    
    for( int i = 0; i < l; i++ ) M.uniformrefinement();

    
    LOG << "Assemble matrices..." << nl;
        
    SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

    MatrixCSR mass = MatrixCSR( massmatrix_scalar );

    LOG << "Prepare data..." << nl;
        
    FloatVector v = mass.createinputvector();
    FloatVector w = mass.createinputvector();

    Float dummy = 0.;

    const int reps = 10;
    const int inner_reps = 100;

    std::clock_t c_start;
    std::clock_t c_end;
    std::clock_t c_sum = 0;

    LOG << "Iterate..." << nl;
        
    for( int i = 0; i < reps; i++ )
    {
        
        v.random();
        w.random();
        
        c_start = std::clock();
        for( int j = 0; j < inner_reps; j++ )
        {
            mass.apply(v,w,1.);
            std::swap(v,w);
        }
        c_end = std::clock();

        c_sum += c_end - c_start; 
        
        Float alpha = v.norm(); 
        v /= alpha;
    }
    
    LOG << "Time passed: " << 1000.0 * (c_sum) / Float(CLOCKS_PER_SEC) << "ms" << nl;
    LOG << "Dummy output to distract optimizer: " << dummy << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
