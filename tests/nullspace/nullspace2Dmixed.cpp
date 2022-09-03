

/**/

#include <fstream>

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
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: Compute a nullspace " << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D Mx = StandardSquare2D();
            
            Mx.check();
            
            Mx.automatic_dirichlet_flags();

            
            
            MeshSimplicial2D M;
            
            for( int i = 0; i < 3; i++ )
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
            

            LOG << "Nullspace computation" << nl;

            ConvergenceTable contable("Mass error");
            
            contable << "#nullvec";
            

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
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                    SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();
                    
                    LOG << "... assemble full matrices" << nl;
            
                    auto physical_mass = volume_incmatrix_t * volume_massmatrix * volume_incmatrix;

                    auto mat_A  = vector_incmatrix_t & vector_massmatrix & vector_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    auto mat_Bt = vector_incmatrix_t & diffmatrix_t & volume_massmatrix & volume_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    
                    auto Z  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                    
                    
                    auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                              + MatrixCSR( vector_incmatrix_t & diffmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix );
                    auto PC = MatrixCSR( volume_incmatrix_t & volume_massmatrix & volume_incmatrix );
                              
                    
                    
//                     auto SystemMatrix = B * inv( A, 1e-10, 0 ) * Bt;
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
                            
                            if( reduced_mass < 1e-6 ) {
                                LOG << "**** The candidate already has very small mass" << nl;
//                                 continue;
                            }
                        }
                        
                        /* reduce the candidate to its nullspace component */
                        {
                            FloatVector rhs( Bt.getdimin(), 0. );
                        
                            FloatVector residual( rhs );
                            
                            for( int t = 0; t < max_number_of_purifications; t++ )
                            {
                                
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
                                
                                if(false)
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
                                
                                assert( candidate.isfinite() );
                                
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
                        
                        if( reduced_mass < 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << nl;
                            continue;
                        }
                        
                        candidate.normalize(mass);
                        
                        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                        
                        LOG << "\t\t\t Numerical residual: " << residual_mass << nl;
                        
                        if( false and residual_mass > 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << nl;
                            continue;
                        }
                        
                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
                    
                        
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
                        
                        // vtk.writeCellVectorData( nullvector, "nullvector H(div)" , 0.1 );
                        
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

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

                contable << nl;
                
                contable.lg();
        
            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}



