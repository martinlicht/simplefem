

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
    
    if(true){

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
                    
        
        LOG << "Nullspace computation" << nl;

        ConvergenceTable contable("Number of nullvectors");
        
        contable << "#nullvec" << nl;
        

        const int min_l = 0; 
        
        const int max_l = 4;
        
        const int min_r = 2; 
        
        const int max_r = 2;
        
        const int max_number_of_candidates = 6;

        const int max_number_of_purifications = 2;

        assert( 0 <= min_l and min_l <= max_l );
        assert( 0 <= min_r and min_r <= max_r );
        
        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ )
        {
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            for( int r = min_r; r <= max_r; r++ )
            {
                
                LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                
                LOG << "...assemble partial matrices" << nl;
        
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
                
                auto Z  = MatrixCSR( B.getdimout(), B.getdimout() ); // zero matrix
                
                auto PA = Conjugation( vector_massmatrix, vector_incmatrix )
                          + 
                          Conjugation( volume_massmatrix, diffmatrix & vector_incmatrix );

                auto PC = Conjugation( volume_massmatrix, volume_incmatrix );

                const auto SystemMatrix = B * inv( A, desired_precision, -1 ) * Bt;
                
                const auto& mass = physical_mass;
                
                
                std::vector<FloatVector> nullvectorgallery;
                
                
                for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
                {
                    FloatVector candidate( Bt.getdimin(), 0. ); 
                    candidate.random(); 
                    candidate.normalize(mass);
                    
                    
                    {
                        for( int s = 0; s < 2; s++ )
                        for( const auto& nullvector : nullvectorgallery ) {
                            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                            candidate = candidate - alpha * nullvector;
                        }
                        
                        Float reduced_mass = candidate.norm(mass);
                        LOG << "\t\t\t Preprocessed mass: " << reduced_mass << nl;
                        
                        if( reduced_mass < mass_threshold_for_small_vectors ) {
                            LOG << "**** The candidate already has very small mass" << nl;
//                                 continue;
                        }
                    }
                    
                    /* reduce the candidate to its nullspace component */
                    {
                        const FloatVector rhs( Bt.getdimin(), 0. );
                    
                        FloatVector residual( rhs );
                        
                        for( int t = 0; t < max_number_of_purifications; t++ )
                        {
                            
                            if( /* DISABLES CODE */ (false) )
                            {

                                const auto PAinv = inv(PA,desired_precision,-1);
                                const auto PCinv = inv(PC,desired_precision,-1);

                                FloatVector  x_A( A.getdimin(),  0. ); 
                                FloatVector& x_C = candidate;
                                
                                const FloatVector  b_A( A.getdimin(),  0. ); 
                                const FloatVector& b_C = rhs; 
                                
                                BlockHerzogSoodhalterMethod( 
                                    x_A, 
                                    x_C, 
                                    b_A, 
                                    b_C, 
                                    -A, Bt, B, Z, 
                                    desired_precision,
                                    -1,
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
                                Z.getA(),   Z.getC(),  Z.getV(), 
                                residual.raw(),
                                desired_precision,
                                -1,
                                desired_precision,
                                -1
                            );
                            
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
                            
                            candidate.normalize( mass );
                            
                            assert( candidate.is_finite() );
                            
                            LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
                            LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << nl;
                            LOG << "\t\t\t (norm eucl) Ax:        " << ( SystemMatrix* candidate ).norm() << nl;
                            LOG << "\t\t\t (norm mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << nl;
                            
                        }
                    }
                    
                    
                    /* Gram-Schmidt */
                    
                    for( int s = 0; s < 2; s++ )
                    for( const auto& nullvector : nullvectorgallery ) {
                        Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                        candidate = candidate - alpha * nullvector;
                    }
                    
                    Float reduced_mass = candidate.norm(mass);
                    LOG << "\t\t\t Reduced mass: " << reduced_mass << nl;
                    
                    if( reduced_mass < mass_threshold_for_small_vectors ) {
                        LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << nl;
                        continue;
                    }
                    
                    candidate.normalize(mass);
                    
                    Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                    
                    LOG << "\t\t\t Numerical residual: " << residual_mass << nl;
                    
                    if( residual_mass > mass_threshold_for_small_vectors ) {
                        LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << nl;
                        continue;
                    }
                    
                    assert( candidate.is_finite() );
                    
                    LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
                
                    
                    nullvectorgallery.push_back( candidate );
                }
                
                
                
                LOG << "How much nullspace are our vectors?" << nl;
                for( const auto& nullvector : nullvectorgallery ) {
                    Float mass_norm = ( SystemMatrix * nullvector ).norm(mass);
                    Assert( mass_norm < mass_threshold_for_small_vectors, mass_norm, mass_threshold_for_small_vectors );
                    // LOGPRINTF( "% 10.5Le\t", (long double)mass_norm );
                    LOG << mass_norm << tab;
                }
                LOG << nl;
                
                LOG << "How orthonormal are our vectors?" << nl;
                for( int n1 = 0; n1 < nullvectorgallery.size(); n1++ ) {
                    for( int n2 = 0; n2 < nullvectorgallery.size(); n2++ ) {
                        auto nullvector1 = nullvectorgallery[n1];
                        auto nullvector2 = nullvectorgallery[n2];
                        Float mass_prod = mass * nullvector1 * nullvector2;
                        // LOGPRINTF( "% 10.5Le\t", (long double)mass_prod );
                        LOG << mass_prod << tab;
                        if( n1 != n2 ) assert( is_numerically_small( mass_prod ) );
                        
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
    
    }
    
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}



