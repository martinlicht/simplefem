

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
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.whitneyincl.hpp"
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

    const int max_number_of_purifications = 2;

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
            
            
            std::vector<FloatVector> nullvectorgallery;
            
            for( int candidate_number = 0; candidate_number < max_number_of_candidates; candidate_number++ )
            {
                
                /* Draw a random vector with unit mass norm as a candidate */

                FloatVector candidate( Bt.getdimin(), 0. ); 
                candidate.random(); 
                candidate.normalize(mass);
                assert( candidate.is_finite() );
                
                /* Orthogonalize that vector against previous nullspace vectors */

                {
                    for( int s = 0; s < 2; s++ )
                    for( const auto& nullvector : nullvectorgallery ) 
                    {
                        Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                        candidate = candidate - alpha * nullvector;
                    }
                }
                
                /* Re-normalize the candidate once again */

                candidate.normalize(mass);

                assert( candidate.is_finite() );
                
                /* reduce the candidate to its nullspace component */
                
                {
                    bool candidate_to_be_discarded = false;
                    
                    for( int t = 0; t < max_number_of_purifications and not candidate_to_be_discarded; t++ )
                    {
                        
                        const FloatVector rhs( Bt.getdimin(), 0. );
                
                        FloatVector residual( rhs );
                        
                        if( /* DISABLES CODE */ (false) )
                        {

                            auto PA = Conjugation( scalar_massmatrix, scalar_incmatrix ); 
                                        + Conjugation( vector_massmatrix, vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix ); 
                            auto PC = Conjugation( vector_massmatrix, vector_incmatrix );
                                        + Conjugation( volume_massmatrix, vector_diffmatrix & vector_incmatrix );
                            
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

                        }
                        
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
                        
                        

                        
                        
                        /*
                        LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm() << nl;
                        LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm( mass ) << nl;
                        LOG << "\t\t\t (eucl) res:       " << residual.norm() << nl;
                        LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << nl;
                        LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
                        LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << nl;
                        LOG << "\t\t\t (eucl) Ax:        " << ( SystemMatrix * candidate ).norm() << nl;
                        LOG << "\t\t\t (mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << nl;
                        LOG << "\t\t\t (eucl) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm() << nl;
                        LOG << "\t\t\t (mass) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( mass ) << nl;
                        */
                        

                        /* Orthogonalize that candidate once again against nullspace vectors */
                        /* Discard if nothing new is found, else re-normalize */
                       
                        for( int s = 0; s < 2; s++ )
                        for( const auto& nullvector : nullvectorgallery ) {
                            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                            candidate = candidate - alpha * nullvector;
                        }
                        
                        Float orthogonalized_candidate_mass = candidate.norm(mass);
                        LOG << "\t\t\t Reduced candidate orthogonalized mass: " << orthogonalized_candidate_mass << nl;

                        if( orthogonalized_candidate_mass < mass_threshold_for_small_vectors ) {
                            // LOGPRINTF("\n\t\t\t !!!!! Discard vector %d/%d because mass is too small!\n\n", candidate_number+1, max_number_of_candidates );
                            candidate_to_be_discarded = true;
                            break;
                        }
                        
                        candidate.normalize(mass);

                        assert( candidate.is_finite() );
                        
                        /*
                        LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
                        LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << nl;
                        */
                        LOG << "\t\t\t (norm eucl) Ax:        " << ( SystemMatrix* candidate ).norm() << nl;
                        LOG << "\t\t\t (norm mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << nl;

                    }

                    // if( candidate_to_be_discarded ) continue;
                }
                
                /* Discard if insufficiently nullspace */ 
                
                Float residual_mass = sqrt( ( C * candidate ).norm_sq( mass ) + ( Bt * candidate ).norm_sq( A ) );
                
                LOG << "\t\t\t Numerical residual: " << residual_mass << nl;
                
                if( residual_mass > mass_threshold_for_small_vectors ) {
                    LOGPRINTF("\n\t\t\t !!!!!Discard vector %d/%d because not nullspace enough!\n\n", candidate_number+1, max_number_of_candidates );
                    continue;
                }
                
                /* Orthogonalize that candidate once again against nullspace vectors */
                
                LOG << "\t\t\t Reduced candidate non-orthogonalized mass: " << candidate.norm(mass) << nl;

                for( int s = 0; s < 2; s++ )
                for( const auto& nullvector : nullvectorgallery ) {
                    Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                    candidate = candidate - alpha * nullvector;
                }
                
                /* Discard if nothing new is found */

                Float orthogonalized_candidate_mass = candidate.norm(mass);
                LOG << "\t\t\t Reduced candidate orthogonalized mass: " << orthogonalized_candidate_mass << nl;

                if( orthogonalized_candidate_mass < mass_threshold_for_small_vectors ) {
                    LOGPRINTF("\n\t\t\t !!!!! Discard vector %d/%d because mass is too small!\n\n", candidate_number+1, max_number_of_candidates );
                    continue;
                }
                
                /* Normalize and accept */

                candidate.normalize(mass);
                
                assert( candidate.is_finite() );
                
                LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
            
                
                nullvectorgallery.push_back( candidate );
            }
            
        
        
            LOG << "How much nullspace are our vectors?" << nl;
            for( const auto& nullvector : nullvectorgallery ) {
                // Float residual_mass = ( SystemMatrix * nullvector ).norm(mass);
                Float residual_mass = sqrt( ( C * nullvector ).norm_sq( mass ) + ( Bt * nullvector ).norm_sq( A ) );
                Assert( residual_mass < mass_threshold_for_small_vectors, residual_mass, mass_threshold_for_small_vectors );
                // LOGPRINTF( "% 10.5Le\t", (long double)residual_mass );
                LOG << residual_mass << nl;
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
            
            
            const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );

            for( const auto& nullvector : nullvectorgallery )
            {
        
                std::fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
    
                VTKWriter vtk( M, fs, getbasename(__FILE__) );

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



