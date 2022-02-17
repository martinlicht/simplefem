

/**/

#include <ostream>
#include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
// #include "../../operators/composed.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/nullspace.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D Mx = StandardSquare2D_tiles3x3();
            
            Mx.automatic_dirichlet_flags();

            // LOG << Mx << nl;

            // Seitenmitten: 2, 11, 19, 31
            Mx.set_flag( 1,  2, SimplexFlagNull );
            Mx.set_flag( 1, 11, SimplexFlagNull );
            Mx.set_flag( 1, 19, SimplexFlagNull );
            Mx.set_flag( 1, 31, SimplexFlagNull );

            // Links: 1, 11, 21 Rechts: 9, 19, 29
            // Mx.set_flag( 1,  1, SimplexFlagNull );
            // Mx.set_flag( 1, 11, SimplexFlagNull );
            // Mx.set_flag( 1, 21, SimplexFlagNull );
            // Mx.set_flag( 1,  9, SimplexFlagNull );
            // Mx.set_flag( 1, 19, SimplexFlagNull );
            // Mx.set_flag( 1, 29, SimplexFlagNull );
            
            // Mx.set_flag( 0, 4, SimplexFlagNull );
            // Mx.set_flag( 0, 8, SimplexFlagNull );
            // Mx.set_flag( 0, 7, SimplexFlagNull );
            // Mx.set_flag( 0,11, SimplexFlagNull );
            



            Mx.check();
            
            
            
            MeshSimplicial2D M;
            
            for( int i = 0; i < 1; i++ )
            {
                auto M2 = Mx;
                M2.getcoordinates().shift( { i * 3.0, 0.0 } );
                M.merge( M2 );
            }
                        
            
            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            

            LOG << "Nullspace computation" << endl;

            ConvergenceTable contable("Mass error");
            
            contable << "#nullvec";

            const Float desired_precision = 100 * machine_epsilon;
            

            const int min_l = 0; 
            
            const int max_l = 4;
            
            const int min_r = 1; 
            
            const int max_r = 1;
            
            const int max_number_of_candidates = 4;

            const int max_number_of_purifications = 1;

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << std::endl;
                    
                    LOG << "...assemble partial matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix vector_elevationmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1);
                    SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r   );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();


                    LOG << "... full matrices" << endl;
            
                    auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                    auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    LOG << vector_elevationmatrix.getdimin() << space << vector_elevationmatrix.getdimout() << nl;

                    auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
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
                    
                    auto SystemMatrix = C + B * inv(A,100*machine_epsilon) * Bt;
                    
                    
                    
                    
                    
                    
                    
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
                                
                                auto& X = SystemMatrix; //B * inv(A,desired_precision) * Bt + C;

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
                                    100*machine_epsilon,
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
                                // Solver.threshold           = desired_precision;
                                // Solver.print_modulo        = 0;
                                // Solver.max_iteration_count = 1 * candidate.getdimension();
                                // Solver.solve( candidate, rhs );
                            
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + X * candidate ).norm() << std::endl;
                                LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + X * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) res:       " << residual.norm() << std::endl;
                                LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) x:         " << candidate.norm() << std::endl;
                                LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) Ax:        " << ( X * candidate ).norm() << std::endl;
                                LOG << "\t\t\t (mass) Ax:        " << ( X * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) b - Ax:    " << ( X * candidate - rhs ).norm() << std::endl;
                                LOG << "\t\t\t (mass) b - Ax:    " << ( X * candidate - rhs ).norm( mass ) << std::endl;
                                
                                candidate.normalize( mass );
                                
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << std::endl;
                                LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t (norm eucl) Ax:        " << ( X * candidate ).norm() << std::endl;
                                LOG << "\t\t\t (norm mass) Ax:        " << ( X * candidate ).norm( mass ) << std::endl;
                                
                                
                                
                            }
                        }
                        
                        
                        /* Is it nullspace at all? */
                        
                        candidate.normalize(mass);
                        
                        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                        
                        LOG << "\t\t\t Numerical residual (after normalizing): " << residual_mass << std::endl;
                        
                        if( residual_mass > 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << std::endl;
                            continue;
                        }
                        
                        /* Gram-Schmidt */
                        
                        for( int s = 0; s < 2; s++ )
                        for( const auto& nullvector : nullvectorgallery ) {
                            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                            candidate = candidate - alpha * nullvector;
                        }
                        
                        Float reduced_mass = candidate.norm(mass);
                        LOG << "\t\t\t Reduced mass: " << reduced_mass << std::endl;
                        LOG << "\t\t\t Numerical residual (after Gram-Schmidt): " << ( SystemMatrix * candidate ).norm(mass) << std::endl;
                        
                        if( reduced_mass < 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << std::endl;
                            continue;
                        }
                        
                        candidate.normalize(mass);

                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector #" << nullvectorgallery.size() + 1 << std::endl;
                    
                        
                        nullvectorgallery.push_back( candidate );
                    }
                    
                    
                    
                    LOG << "How much nullspace are our vectors?" << nl;
                    for( const auto& nullvector : nullvectorgallery ) {
                        LOGPRINTF( "% 10.5e\t", ( SystemMatrix * nullvector ).norm(mass) );
                        // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << ( SystemMatrix * nullvector ).norm(mass) << tab;
                    }
                    LOG << nl;
                    
                    LOG << "How orthonormal are our vectors?" << nl;
                    for( const auto& nullvector1 : nullvectorgallery ) {
                        for( const auto& nullvector2 : nullvectorgallery ) {
                            LOGPRINTF( "% 10.5e\t", mass * nullvector1 * nullvector2 );
                            // LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << mass * nullvector1 * nullvector2 << tab;
                        }
                        LOG << nl;
                    }
                    
                    
                    
                    contable << static_cast<Float>(nullvectorgallery.size());  

                    if( r == 1 )
                    for( const auto& nullvector : nullvectorgallery )
                    {
                
                        fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
            
                        VTKWriter vtk( M, fs, getbasename(__FILE__) );
                        vtk.writeCoordinateBlock();
                        vtk.writeTopDimensionalCells();
                        
                        // vtk.writeCellVectorData( nullvector, "nullvector H(curl)" , 0.1 );
                        
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
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << std::endl;
//                         
//                         }
//                         
//                         
//                         
//                         contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
//                         
//                         
//                         if( r == 1 ) {
//                     
//                             fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
//                 
//                             VTKWriter vtk( M, fs, getbasename(__FILE__) );
//                             vtk.writeCoordinateBlock();
//                             vtk.writeTopDimensionalCells();
//                             
//                             vtk.writeVertexScalarData( sol,  "data1" , 1.0 );
// //                             vtk.writeVertexScalarData( sol2, "data2" , 1.0 );
//                             // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
//                             
//                             fs.close();
//                     
//                         }
// 
//                             
//                             
//                         contable << nl;
//                         
//                         contable.lg( false );
// 
//                     }
                    
                }

                LOG << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();

                contable << nl;
                
                contable.lg();
        
            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}



