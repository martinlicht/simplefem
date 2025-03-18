

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
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"


// using namespace std;

const Float mass_threshold_for_small_vectors = 1e-6;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Compute a nullspace " << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D Mx = StandardSquare2D();
    
    Mx.check();
    
    Mx.automatic_dirichlet_flags();

    
    
    MeshSimplicial2D M;
    
    for( int i = 0; i < 3; i++ )
    {
        auto M2 = Mx;
        M2.getCoordinates().shift( FloatVector{ i * 3.0, 0.0 } );
        M.merge( M2 );
    }
                
    
    const int min_l = 0; 
    
    const int max_l = 4;
    
    const int min_r = 2; 
    
    const int max_r = 2;
    
    const int max_number_of_candidates = 6;

    // const int max_number_of_purifications = 2;

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    
    ConvergenceTable contable("Nullvectors found");
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
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        for( int r = min_r; r <= max_r; r++ )
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
            
            LOG << "... assemble partial matrices" << nl;
    
            auto vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   ) );
            
            auto volume_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ) );

            auto diffmatrix   = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r ) );
            auto diffmatrix_t = diffmatrix.getTranspose();

            auto vector_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   ) );
            auto vector_incmatrix_t = vector_incmatrix.getTranspose();

            auto volume_incmatrix   = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 ) );
            auto volume_incmatrix_t = volume_incmatrix.getTranspose();
            
            LOG << "... assemble full matrices" << nl;
    
            auto physical_mass = Conjugation( volume_massmatrix, volume_incmatrix );

            auto A  = Conjugation( vector_massmatrix, vector_incmatrix );
            
            auto Bt = vector_incmatrix_t & diffmatrix_t & volume_massmatrix & volume_incmatrix; // upper right
            
            auto B  = Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower left
            
            auto C  = MatrixCSR( B.getdimout(), B.getdimout() ); // zero matrix
            
            auto PA = Conjugation( vector_massmatrix, vector_incmatrix )
                        + 
                        Conjugation( volume_massmatrix, diffmatrix & vector_incmatrix );

            auto PC = Conjugation( volume_massmatrix, volume_incmatrix );

            const auto SystemMatrix = B * inv( A, desired_precision, -1 ) * Bt;
            
            const auto& mass = physical_mass;

            
            
            
            auto purifier = [&]( FloatVector& candidate ){
                
                FloatVector rhs( candidate.getdimension(), 0. );
                
                // FloatVector residual( rhs );
                // HodgeConjugateResidualSolverCSR_SSOR(
                //     B.getdimout(), 
                //     A.getdimout(), 
                //     candidate.raw(), 
                //     rhs.raw(), 
                //     A.getA(),   A.getC(),  A.getV(), 
                //     B.getA(),   B.getC(),  B.getV(), 
                //     Bt.getA(), Bt.getC(), Bt.getV(), 
                //     C.getA(),   C.getC(),  C.getV(), 
                //     residual.raw(),
                //     desired_precision,
                //     -3,
                //     desired_precision,
                //     -3
                // );

                const auto PAinv = inv(PA,desired_precision,-3);
                const auto PCinv = inv(PC,desired_precision,-3);
                FloatVector  x_A( A.getdimin(),  0. ); 
                FloatVector& x_C = candidate;
                const FloatVector  b_A( A.getdimin(),  0. ); 
                const FloatVector& b_C = rhs; 
                BlockHerzogSoodhalterMethod( 
                    x_A, 
                    x_C, 
                    b_A, 
                    b_C, 
                    -A, Bt, B, C, 
                    desired_precision,
                    -3,
                    PAinv, PCinv
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

            // HodgeConjugateResidualSolverCSR_SSOR(
            //     B.getdimout(), 
            //     A.getdimout(), 
            //     candidate.raw(), 
            //     rhs.raw(), 
            //     A.getA(),   A.getC(),  A.getV(), 
            //     B.getA(),   B.getC(),  B.getV(), 
            //     Bt.getA(), Bt.getC(), Bt.getV(), 
            //     C.getA(),   C.getC(),  C.getV(), 
            //     residual.raw(),
            //     desired_precision,
            //     -3,
            //     desired_precision,
            //     -3
            // );

            // ConjugateResidualSolverCSR( 
            //     candidate.getdimension(), 
            //     candidate.raw(), 
            //     rhs.raw(), 
            //     C.getA(), C.getC(), C.getV(),
            //     residual.raw(),
            //     desired_precision,
            //     0
            // );

            // {
            //     const FloatVector newrhs = SystemMatrix * candidate;
            //     FloatVector input = candidate; input.zero();
            //     ConjugateResidualMethod solver( SystemMatrix );
            //     solver.precision           = desired_precision;
            //     solver.print_modulo        = -2;
            //     solver.max_iteration_count = 1 * candidate.getdimension();
            //     solver.solve( input, newrhs );
            //     candidate = candidate - input;
            //     assert( candidate.is_finite() );
            // }

            // {
            //     auto PA = Conjugation( scalar_massmatrix, scalar_incmatrix ); 
            //                 + Conjugation( vector_massmatrix, vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix ); 
            //     auto PC = Conjugation( vector_massmatrix, vector_incmatrix );
            //                 + Conjugation( volume_massmatrix, vector_diffmatrix & vector_incmatrix );
            //     const auto PAinv = inv(PA,desired_precision,-3);
            //     const auto PCinv = inv(PC,desired_precision,-3);
            //     FloatVector x_A( A.getdimin(),  0. ); 
            //     FloatVector x_C( C.getdimin(),  0. );
            //     const FloatVector b_A = -A * x_A + Bt * candidate;
            //     const FloatVector b_C =  B * x_A +  C * candidate;
            //     for( int i = 0; i < 3; i++ )
            //     BlockHerzogSoodhalterMethod( 
            //         x_A, 
            //         x_C, 
            //         b_A, 
            //         b_C, 
            //         -A, Bt, B, C, 
            //         desired_precision,
            //         0,
            //         PAinv, PCinv
            //     );
            //     candidate = candidate - x_C;
            // }
            
            // {
            //     const auto& X = Block2x2Operator( -A, Bt, B, C );
            //     FloatVector zeroA = A.createinputvector(); zeroA.zero();
            //     FloatVector foo = candidate; foo.random();
            //     FloatVector input = concatFloatVector( zeroA, foo );
            //     const FloatVector newrhs = X * input;
            //     LOG << input.norm() << space << newrhs.norm() << nl;
            //     FloatVector newinput = input;
            //     MinimumResidualMethod solver( X );
            //     solver.precision           = desired_precision;
            //     solver.print_modulo        = 1;
            //     solver.max_iteration_count = input.getdimension();
            //     solver.solve( newinput, newrhs );
            // }
            
            
            
        
        
            
            
            
            
            LOG << "How much nullspace are our vectors?" << nl;
            for( const auto& nullvector : nullvectorgallery ) {
                // Float residual_mass = ( SystemMatrix * nullvector ).norm(mass);
                Float residual_mass = sqrt( ( C * nullvector ).norm_sq( mass ) + ( Bt * nullvector ).norm_sq( A ) );
                Assert( residual_mass < mass_threshold_for_small_vectors, residual_mass, mass_threshold_for_small_vectors );
                LOGPRINTF( "% 10.5le\t", (double)(safedouble)residual_mass );
            }
            LOG << nl;
            
            LOG << "How orthonormal are our vectors?" << nl;
            for( int n1 = 0; n1 < nullvectorgallery.size(); n1++ ) {
                for( int n2 = 0; n2 < nullvectorgallery.size(); n2++ ) {
                    auto nullvector1 = nullvectorgallery[n1];
                    auto nullvector2 = nullvectorgallery[n2];
                    Float mass_prod = mass * nullvector1 * nullvector2;
                    LOGPRINTF( "% 10.5le\t", (double)(safedouble)mass_prod );
                    // LOG << mass_prod << tab;
                    if( n1 != n2 ) 
                        assert( is_numerically_small( mass_prod ) );
                    else
                        assert( is_numerically_one( mass_prod ) );
                    
                }
                LOG << nl;
            }
            
            
            contable << static_cast<Float>(nullvectorgallery.size());
            
            
            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r-1 );

            for( const auto& nullvector : nullvectorgallery )
            {
        
                std::fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
    
                VTKWriter vtk( M, fs, getbasename(__FILE__) );
                
                auto reduced_nullvector = interpol_matrix * volume_incmatrix * nullvector;

                vtk.write_cell_scalar_data_barycentricvolumes( reduced_nullvector, "nullvector_L2" , 1.0 );
                
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



