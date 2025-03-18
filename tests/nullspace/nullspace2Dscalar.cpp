

/**/

#include <fstream>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/nullspace.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.interpol.hpp"


// using namespace std;

const Float mass_threshold_for_small_vectors = 1e-6;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Nullspace computation (2D) scalar" << nl;

    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D Mx = StandardSquare2D();
    
    MeshSimplicial2D M;
    
    for( int i = 0; i < 3; i++ )
    {
        auto M2 = Mx;
        M2.getCoordinates().shift( FloatVector{ i * 3.0, 0.0 } );
        M.merge( M2 );
    }
                
    M.check();
    
    
    bool do_sullivan = true;
    
    bool do_whitney  = true;
    
    const int min_l = 0; 
    
    const int max_l = 6;
    
    const int min_r = 1; 
    
    const int max_r = 2;
    
    const int max_number_of_candidates = 6;

    const int max_number_of_purifications = 2;

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    
    ConvergenceTable contable("Nullvectors found");
    contable.display_convergence_rates = false;

    for( int r = min_r; r <= max_r; r++ )
    for( int b = 0; b <= 1; b++ )
    {
        if( b == 0 and not do_sullivan ) continue;
        if( b == 1 and not do_whitney  ) continue;
        contable << printf_into_string("#nullvec%c%i", b?'W':'S', r );
    }
    contable << nl;
    
    
    LOG << "Nullspace computation" << nl;
    
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ )
    {
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        for( int r = min_r; r <= max_r; r++ )
        for( int b = 0; b <= 1; b++ )
        {
            
            if( b == 0 and not do_sullivan ) continue;
            if( b == 1 and not do_whitney  ) continue;
            
            LOG << "Polynomial degree: " << r << "/" << max_r << " using " << (b==0?"Sullivan":"Whitney") << " forms" << nl;
            
            LOG << "... assemble matrices" << nl;
    
            const auto scalar_massmatrix = MatrixCSR(FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r ));
            
            LOG << "... assemble vector mass matrix" << nl;
    
            const auto vector_massmatrix = MatrixCSR(FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 ));
            
            LOG << "... assemble differential matrix and transpose" << nl;

            const auto diffmatrix = MatrixCSR(FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r ));

            LOG << "... assemble inclusion matrix and transpose" << nl;
    
            const auto incmatrix = ( b == 0 ) ? 
                                   MatrixCSR(FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r )) 
                                   : 
                                   MatrixCSR(FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r ));

            LOG << "... assemble stiffness matrix" << nl;
            
            const auto stiffness = Conjugation( (vector_massmatrix), (diffmatrix) & (incmatrix) );
            
            const auto physical_mass = Conjugation( (scalar_massmatrix), (incmatrix) );
            
            const auto& SystemMatrix = stiffness;
            
            const auto& mass = physical_mass;
            
            const auto precon = InverseDiagonalPreconditioner(SystemMatrix);

            auto purifier_cpp = [&]( FloatVector& candidate ){
                ConjugateResidualMethod CRM(SystemMatrix);
                CRM.verbosity = IterativeSolver::VerbosityLevel::startandfinish;

                FloatVector solution( candidate.getdimension(), 0. );
                CRM.solve( solution, SystemMatrix * candidate );
                candidate -= solution;
            };
            
            auto purifier_csr = [&]( FloatVector& candidate ){
                
                FloatVector solution( candidate.getdimension(), 0. );
                
                const FloatVector rhs( SystemMatrix.getdimin(), 0. );
                FloatVector residual( rhs );
                     
                ConjugateResidualSolverCSR( 
                    candidate.getdimension(), 
                    candidate.raw(), 
                    rhs.raw(), 
                    SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                    residual.raw(),
                    desired_precision,
                    -3
                );
                
                // candidate -= solution;
            };
            
            auto purifier_csr2 = [&]( FloatVector& candidate ){
                
                FloatVector solution( candidate.getdimension(), 0. );
                
                const FloatVector rhs = SystemMatrix * candidate;
                FloatVector residual( rhs );
                     
                ConjugateResidualSolverCSR( 
                    solution.getdimension(), 
                    solution.raw(), 
                    rhs.raw(), 
                    SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                    residual.raw(),
                    desired_precision,
                    -3
                );
                
                candidate -= solution;
            };
            
            std::vector<FloatVector> nullvectorgallery = computeNullspace(
                SystemMatrix,
                physical_mass,
                physical_mass,
                max_number_of_candidates,
                //
                mass_threshold_for_small_vectors,
                mass_threshold_for_small_vectors,
                purifier_csr2
            );

            
            
            
            contable << static_cast<Float>(nullvectorgallery.size());

            
            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );

            for( const auto& nullvector : nullvectorgallery )
            {
        
                std::fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
    
                VTKWriter vtk( M, fs, getbasename(__FILE__) );
                
                auto reduced_nullvector = interpol_matrix * incmatrix * nullvector;

                vtk.write_cell_scalar_data( reduced_nullvector,  "nullvector" , 1.0 );
                
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


