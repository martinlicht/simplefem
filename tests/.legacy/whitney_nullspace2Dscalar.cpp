

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
#include "../../solver/sparsesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.interpol.hpp"


// using namespace std;

const Float mass_threshold_for_small_vectors = 1e-6;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: Nullspace computation (2D) scalar" << nl;
    
    if(true){

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
        
        
        LOG << "Nullspace computation" << nl;

        
        ConvergenceTable contable("Nullvectors found");
        
        contable << "#nullvec" << nl;
                    

        const int min_l = 0; 
        
        const int max_l = 6;
        
        const int min_r = 1; 
        
        const int max_r = 1;
        
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
                
                LOG << "...assemble matrices" << nl;
        
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                LOG << "...assemble vector mass matrix" << nl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                
                LOG << "...assemble differential matrix and transpose" << nl;

                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                LOG << "...assemble inclusion matrix and transpose" << nl;
        
                SparseMatrix incmatrix = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );

                SparseMatrix incmatrix_t = incmatrix.getTranspose();

                LOG << "...assemble stiffness matrix" << nl;

                auto stiffness = Conjugation( MatrixCSR(vector_massmatrix), MatrixCSR(diffmatrix) & MatrixCSR(incmatrix) );
                
                auto physical_mass = Conjugation( MatrixCSR(scalar_massmatrix), MatrixCSR(incmatrix) );
                
                auto physical_mass = MatrixCSR( incmatrix_t & ( scalar_massmatrix & incmatrix ) );
                
                const auto& SystemMatrix = stiffness;
                
                std::vector<FloatVector> nullvectorgallery;
                
                const auto& mass = physical_mass;
                
                for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
                {
                    
                    FloatVector candidate( SystemMatrix.getdimin(), 0. ); 
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
                        for( int t = 0; t < max_number_of_purifications; t++ )
                        {
                            
                            const FloatVector rhs( SystemMatrix.getdimin(), 0. );
                            FloatVector residual( rhs );
                        
                            ConjugateResidualSolverCSR( 
                                candidate.getdimension(), 
                                candidate.raw(), 
                                rhs.raw(), 
                                SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            
//                                 FloatVector zero( candidate.getdimension(), 0. );
//                                 FloatVector residual( candidate.getdimension(), 0. );
//                                 
//                                 ConjugateResidualSolverCSR( 
//                                     zero.getdimension(), 
//                                     zero.raw(), 
//                                     candidate.raw(), 
//                                     SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                     residual.raw(),
//                                     desired_precision,
//                                     0
//                                 );
//                                 candidate = residual;

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
    
    }
    
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}


