

/**/

#include <fstream>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/nullspace.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.inclwhitney.hpp"
#include "../../fem/global.interpol.hpp"


// using namespace std;

const Float mass_threshold_for_small_vectors = 1e-6;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Nullspace computation (2D) Hodge-Laplacian" << nl;

    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D Mx = StandardSquare2D_tiles3x3();
    
    Mx.automatic_dirichlet_flags();

    // LOG << Mx << nl;

    // Seitenmitten: 2, 11, 19, 31
    Mx.set_flag( 1,  2, SimplexFlag::SimplexFlagNull );
    Mx.set_flag( 1, 11, SimplexFlag::SimplexFlagNull );
    Mx.set_flag( 1, 19, SimplexFlag::SimplexFlagNull );
    Mx.set_flag( 1, 31, SimplexFlag::SimplexFlagNull );

    // Links: 1, 11, 21 Rechts: 9, 19, 29
    // Mx.set_flag( 1,  1, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 1, 11, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 1, 21, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 1,  9, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 1, 19, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 1, 29, SimplexFlag::SimplexFlagNull );
    
    // Mx.set_flag( 0, 4, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 0, 8, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 0, 7, SimplexFlag::SimplexFlagNull );
    // Mx.set_flag( 0,11, SimplexFlag::SimplexFlagNull );
    
    Mx.check();
    
    
    
    MeshSimplicial2D M;
    
    for( int i = 0; i < 1; i++ )
    {
        auto M2 = Mx;
        M2.getCoordinates().shift( FloatVector{ i * 3.0, 0.0 } );
        M.merge( M2 );
    }
                
    
    
    const Float desired_precision = 100 * machine_epsilon;
    

    const int min_l = 0; 
    
    const int max_l = 4;
    
    const int min_r = 1; 
    
    const int max_r = 2;
    
    const int max_number_of_candidates = 4;

    // const int max_number_of_purifications = 2;

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    
    ConvergenceTable contable("Nullvectors found");
    contable.display_convergence_rates = false;
    
    for( int r = min_r; r <= max_r; r++ )
    {
        contable << printf_into_string("#nullvec%i", r );
    }
    contable << nl;
    
    
    LOG << "Nullspace computation" << nl;

    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ )
    {
        
        LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        for( int r = min_r; r <= max_r; r++ )
        {
            
            LOG << "Level: " << min_l << " <= " << l << " <= " << max_l << nl;
            LOG << "Polynomial degree: " <<  min_r << " <= " << r << " <= " << max_r << nl;
            
            LOG << "... assemble partial matrices" << nl;
    
            const auto scalar_massmatrix = MatrixCSR(FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   ));
            const auto vector_massmatrix = MatrixCSR(FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   ));
            const auto volume_massmatrix = MatrixCSR(FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ));

            const auto vector_elevationmatrix   = MatrixCSR(FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1));
            const auto vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

            const auto scalar_incmatrix   = MatrixCSR(FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r   ));
            const auto scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            const auto vector_incmatrix   = MatrixCSR(FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   ));
            const auto vector_incmatrix_t = vector_incmatrix.getTranspose();

            const auto scalar_diffmatrix   = MatrixCSR(FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r ));
            const auto scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            const auto vector_diffmatrix   = MatrixCSR(FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r ));
            const auto vector_diffmatrix_t = vector_diffmatrix.getTranspose();


            LOG << "... full matrices" << nl;
    
            auto mass = Conjugation( vector_massmatrix, vector_incmatrix );

            const auto A  = Conjugation( scalar_massmatrix, scalar_incmatrix );
            
            const auto Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
            
            const auto B = Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
            
            const auto C = Conjugation( volume_massmatrix, vector_diffmatrix & vector_incmatrix );
            
            auto SystemMatrix = C + B * inv(A,100*machine_epsilon,-3) * Bt;
            
            
            auto purifier = [&]( FloatVector& candidate ){
                
                FloatVector rhs( candidate.getdimension(), 0. );
                
                FloatVector residual( rhs );
                     
                HodgeConjugateResidualSolverCSR_SSOR(
                    B.getdimout(), 
                    A.getdimout(), 
                    candidate.raw(), 
                    rhs.raw(), 
                    A.getA(),   A.getC(),  A.getV(), 
                    B.getA(),   B.getC(),  B.getV(), 
                    Bt.getA(), Bt.getC(), Bt.getV(), 
                    C.getA(),   C.getC(),  C.getV(), 
                    residual.raw(),
                    desired_precision,
                    -3,
                    desired_precision,
                    -3
                );

            };
            
            std::vector<FloatVector> nullvectorgallery = computeNullspace(
                Block2x2Operator( 
                    ZeroOperator(A.getdimout(),0), Bt, 
                    ZeroOperator(B.getdimout(),0), C   ),
                mass,
                Block2x2Operator( 
                    A,                                        ZeroOperator(Bt.getdimout(),Bt.getdimin()),
                    ZeroOperator(B.getdimout(),B.getdimin()), mass                                        ),
                max_number_of_candidates,
                //
                mass_threshold_for_small_vectors,
                mass_threshold_for_small_vectors,
                purifier
            );

            
            contable << static_cast<Float>(nullvectorgallery.size());
            
            
            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );

            for( const auto& nullvector : nullvectorgallery )
            {
        
                std::fstream fs( get_available_filename(get_basename(__FILE__)), std::fstream::out );
    
                VTKWriter vtk( M, fs, get_basename(__FILE__) );

                auto reduced_nullvector = interpol_matrix * vector_incmatrix * nullvector;

                vtk.write_cell_vector_data_barycentricgradients( reduced_nullvector, "nullvector_Hcurl" , 1.0 );
                
                fs.close();
        
            } 

        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

        contable << nl;
        
        contable.lg();

    } 
    
    contable.lg();
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}



