

/**/

#include <ostream>
#include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D Mx = StandardSquare2D();
            
            Mx.check();
            
            Mx.automatic_dirichlet_flags();

            
            
            MeshSimplicial2D M;
            
            for( int i = 0; i < 1; i++ )
            {
                auto M2 = Mx;
                M2.getcoordinates().shift( FloatVector{ i * 3.0, 0.0 } );
                M.merge( M2 );
            }
                        
            
            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            

            LOG << "Nullspace computation" << nl;

            ConvergenceTable contable("Mass error");
            
            contable << "#nullvec";
            

            const int min_l = 0; 
            
            const int max_l = 3;
            
            const int min_r = 2; 
            
            const int max_r = 2;
            
            const int max_number_of_candidates = 4;

            const int max_number_of_purifications = 1;

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
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                    SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();
                    
                    LOG << "... assemble full matrices" << nl;
            
                    auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                    auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                    mat_C.sortandcompressentries();
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    auto C  = MatrixCSR( mat_C  );
                    
                    auto Z  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                    
                    auto SystemMatrix = C + B * inv(A,desired_precision, 1) * Bt;
                    
                    
                    
                    
                    
                    
                    
                    std::vector<FloatVector> nullvectorgallery;
                    
                    
                    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
                    {
                        
                        FloatVector candidate( Bt.getdimin(), 0. ); 
                        candidate.random(); 
                        candidate.normalize(mass);
                        
                        /* reduce the candidate to its nullspace component */
                        {
                            FloatVector rhs( Bt.getdimin(), 0. );
                        
                            FloatVector residual( rhs );
                            
                            for( int t = 0; t < max_number_of_purifications; t++ )
                            {
                                
                                auto X = B * inv(A,desired_precision) * Bt + C;

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
                                    0,
                                    desired_precision,
                                    -1
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
                                
                                // ConjugateResidualMethod Solver( X );
                                // Solver.tolerance           = desired_precision;
                                // Solver.print_modulo        = 100;
                                // Solver.max_iteration_count = 4 * candidate.getdimension();
                                // Solver.solve( candidate, rhs );
                            
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + X * candidate ).norm() << nl;
                                LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + X * candidate ).norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) res:       " << residual.norm() << nl;
                                LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
                                LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) Ax:        " << ( X * candidate ).norm() << nl;
                                LOG << "\t\t\t (mass) Ax:        " << ( X * candidate ).norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) b - Ax:    " << ( X * candidate - rhs ).norm() << nl;
                                LOG << "\t\t\t (mass) b - Ax:    " << ( X * candidate - rhs ).norm( mass ) << nl;
                                
                                candidate.normalize( mass );
                                
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
                                LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << nl;
                                LOG << "\t\t\t (norm eucl) Ax:        " << ( X * candidate ).norm() << nl;
                                LOG << "\t\t\t (norm mass) Ax:        " << ( X * candidate ).norm( mass ) << nl;
                                
                                
                                
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
                        
                        if( reduced_mass < 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << nl;
                            continue;
                        }
                        
                        candidate.normalize(mass);
                        
                        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                        
                        LOG << "\t\t\t Numerical residual: " << residual_mass << nl;
                        
                        if( residual_mass > 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << nl;
                            continue;
                        }
                        
                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector #" << nullvectorgallery.size() + 1 << nl;
                    
                        
                        nullvectorgallery.push_back( candidate );
                    }
                    
                    
                    
                    LOG << "How much nullspace are our vectors?" << nl;
                    for( const auto& nullvector : nullvectorgallery ) {
                        LOGPRINTF( "% 10.5e\t", ( SystemMatrix * nullvector ).norm(mass) );
                        // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << ( SystemMatrix * nullvector ).norm(mass) << tab;
                    }
                    LOG << nl;
                    
                    LOG << "How orthonormal are our vectors?" << nl;
                    for( int n1 = 0; n1 < nullvectorgallery.size(); n1++ ) {
                        for( int n2 = 0; n2 < nullvectorgallery.size(); n2++ ) {
                            auto nullvector1 = nullvectorgallery[n1];
                            auto nullvector2 = nullvectorgallery[n2];
                            Float mass_prod = mass * nullvector1 * nullvector2;
                            LOGPRINTF( "% 10.5e\t", mass_prod );
                            if( n1 != n2 ) assert( is_numerically_small( mass_prod ) );
                            
                        }
                        LOG << nl;
                    }
                    
                    
                    
                    contable << static_cast<Float>(nullvectorgallery.size());   
                    
                    
                    const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );

                    for( const auto& nullvector : nullvectorgallery )
                    {
                
                        fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
            
                        VTKWriter vtk( M, fs, getbasename(__FILE__) );
                        // vtk.writeCoordinateBlock();
                        // vtk.writeTopDimensionalCells();

                        auto reduced_nullvector = interpol_matrix * vector_incmatrix * nullvector;

                        vtk.writeCellVectorData_barycentricgradients( reduced_nullvector, "nullvector_Hcurl" , 0.1 );
                        
                        fs.close();
                
                    } 
                    
//                     {
// 
//                         FloatVector sol( opr.getdimin(), 0. ); sol.random(); sol.normalize(mass);
//                         
//                         assert( sol.isfinite() );
//                         
//                         FloatVector rhs( opr.getdimin(), 0. );
//                         
//                         FloatVector residual( rhs );
//                         
//                         for( int t = 0; t < 3; t++ ) {
//                             
//                             ConjugateResidualSolverCSR( 
//                                 sol.getdimension(), 
//                                 sol.raw(), 
//                                 rhs.raw(), 
//                                 mat.getA(), mat.getC(), mat.getV(),
//                                 residual.raw(),
//                                 desired_precision,
//                                 1
//                             );
//                             sol.normalize( mass );
//                             
//                             assert( sol.isfinite() );
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << nl;
//                             LOG << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << nl;
//                             LOG << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << nl;
//                         
//                         }
//                         
//                         
//                         
//                         contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
//                         
//                         
// 
//                             
//                             
//                         contable << nl;
//                         
//                         contable.lg( false );
// 
//                     }
                    
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



