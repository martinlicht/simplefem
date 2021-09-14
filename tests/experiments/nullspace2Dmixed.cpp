

/**/

#include <iostream>
#include <fstream>
#include <iomanip>

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
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
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

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Compute a nullspace" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
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
            

            LOG << "Nullspace computation" << endl;

            ConvergenceTable contable;
            
            contable << "#nullvec";
            

            int min_l = 0; 
            
            int max_l = 4;
            
            int min_r = 2; 
            
            int max_r = 2;
            
            int max_number_of_candidates = 6;

            int max_number_of_purifications = 1;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << std::endl;
                    
                    LOG << "...assemble matrices" << endl;
            
                    LOG << "... assemble matrices" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix vector_massmatrix_inv = FEECBrokenMassMatrix_cellwiseinverse( M, M.getinnerdimension(), 1, r   );
                    
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix volume_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r-1 );
                    SparseMatrix volume_incmatrix_t = volume_incmatrix.getTranspose();
                    
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
                    
                    
//                     auto SystemMatrix = B * inv( A, 1e-10, 0 ) * Bt;
                    const auto SystemMatrix = B * inv( A, 1000*machine_epsilon, -1 ) * Bt;
                    
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
                            LOG << "\t\t\t Preprocessed mass: " << reduced_mass << std::endl;
                            
                            if( reduced_mass < 1e-6 ) {
                                LOG << "**** The candidate already has very small mass" << std::endl;
//                                 continue;
                            }
                        }
                        
                        /* reduce the candidate to its nullspace component */
                        {
                            FloatVector rhs( Bt.getdimin(), 0. );
                        
                            FloatVector residual( rhs );
                            
                            for( int t = 0; t < max_number_of_purifications; t++ )
                            {
                                
                                HodgeConjugateResidualSolverCSR_SSOR(
                                    B.getdimout(), 
                                    A.getdimout(), 
                                    candidate.raw(), 
                                    rhs.raw(), 
                                    A.getA(),   A.getC(),  A.getV(), 
                                    B.getA(),   B.getC(),  B.getV(), 
                                    Bt.getA(), Bt.getC(), Bt.getV(), 
                                    residual.raw(),
                                    100 * machine_epsilon,
                                    -1,
                                    desired_precision,
                                    -1
                                );
                                
                                LOG << "\t\t\t x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( mass ) << std::endl;
                                
                                candidate.normalize( mass );
                                
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( mass ) << std::endl;
                                
                            }
                        }
                        
                        
                        /* Gram-Schmidt */
                        
                        for( int s = 0; s < 2; s++ )
                        for( const auto& nullvector : nullvectorgallery ) {
                            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                            candidate = candidate - alpha * nullvector;
                        }
                        
                        Float reduced_mass = candidate.norm(mass);
                        LOG << "\t\t\t Reduced mass: " << reduced_mass << std::endl;
                        
                        if( reduced_mass < 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << std::endl;
                            continue;
                        }
                        
                        candidate.normalize(mass);
                        
                        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                        
                        LOG << "\t\t\t Numerical residual: " << residual_mass << std::endl;
                        
                        if( false and residual_mass > 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << std::endl;
                            continue;
                        }
                        
                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << std::endl;
                    
                        
                        nullvectorgallery.push_back( candidate );
                    }
                    
                    
                    
                    LOG << "How much nullspace are our vectors?" << nl;
                    for( const auto& nullvector : nullvectorgallery ) {
                        LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << ( SystemMatrix * nullvector ).norm(mass) << tab;
                    }
                    LOG << nl;
                    
                    LOG << "How orthonormal are our vectors?" << nl;
                    for( const auto& nullvector1 : nullvectorgallery ) {
                        for( const auto& nullvector2 : nullvectorgallery ) {
                            LOG << std::showpos << std::scientific << std::setprecision(5) << std::setw(10) << mass * nullvector1 * nullvector2 << tab;
                        }
                        LOG << nl;
                    }
                    
                    
                    
                    contable << static_cast<Float>(nullvectorgallery.size());   
                    
                    
                    
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
//                                 100 * machine_epsilon,
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
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}




//                      {
// 
//                         FloatVector sol_original( opr.getdimin(), 0. );
//                         FloatVector rhs_original( opr.getdimin(), 0. );
//                         
//                         sol_original.random(); sol_original.normalize( mass );
//                         rhs_original.random(); rhs_original.normalize( mass );
//                         
//                         if(false)
//                         {
//                             LOG << "Filter out from x (CGM)" << endl;
//                         
//                             FloatVector sol( sol_original );
//                             FloatVector rhs( rhs_original.getdimension(), 0. );
//                             FloatVector residual( rhs );
//                             
//                             assert( sol.isfinite() );
//                             
//                             
//                             ConjugateGradientSolverCSR( 
//                                 sol.getdimension(), 
//                                 sol.raw(), 
//                                 rhs.raw(), 
//                                 mat.getA(), mat.getC(), mat.getV(),
//                                 residual.raw(),
//                                 machine_epsilon,
//                                 1
//                             );
//                             sol.normalize( mass );
//                             
//                             assert( sol.isfinite() );
//                             
//                             LOG << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
//                             
//                             
//                             FloatVector sol2( sol_original );
//                             sol2.random();
//                             FloatVector rhs2( rhs_original.getdimension(), 0. );
//                             FloatVector residual2( rhs );
//                             
//                             ConjugateGradientSolverCSR( 
//                                 sol2.getdimension(), 
//                                 sol2.raw(), 
//                                 rhs2.raw(), 
//                                 mat.getA(), mat.getC(), mat.getV(),
//                                 residual2.raw(),
//                                 machine_epsilon,
//                                 1
//                             );
//                             sol2.normalize( mass );
//                             
//                             
//                             
//                             
//                             if( r == 1 ) {
//                         
//                                 fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
//                     
//                                 VTKWriter vtk( M, fs, getbasename(__FILE__) );
//                                 vtk.writeCoordinateBlock();
//                                 vtk.writeTopDimensionalCells();
//                                 
//                                 vtk.writeVertexScalarData( sol,  "data1" , 1.0 );
//                                 vtk.writeVertexScalarData( sol2, "data2" , 1.0 );
//                                 // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
//                                 
//                                 fs.close();
//                         
//                             }
// 
//                             
//                             
//                         }
// 
//                         if(false)
//                         {
//                             LOG << "Filter out from x (CRM)" << endl;
//                         
//                             FloatVector sol( sol_original );
//                             FloatVector rhs( rhs_original.getdimension(), 0. );
//                             FloatVector residual( rhs );
//                             
//                             ConjugateResidualSolverCSR( 
//                                 sol.getdimension(), 
//                                 sol.raw(), 
//                                 rhs.raw(), 
//                                 mat.getA(), mat.getC(), mat.getV(),
//                                 residual.raw(),
//                                 machine_epsilon,
//                                 0
//                             );
//                             sol.normalize( mass );
//                             
//                             LOG << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
//                         }
// 
//                         if(false)
//                         {
//                             LOG << "Filter out from b" << endl;
//                         
//                             FloatVector sol( sol_original.getdimension(), 0. );
//                             FloatVector rhs( rhs_original );
//                             FloatVector residual( rhs );
//                             
//                             ConjugateResidualSolverCSR_textbook( 
//                                 sol.getdimension(), 
//                                 sol.raw(), 
//                                 rhs.raw(), 
//                                 mat.getA(), mat.getC(), mat.getV(),
//                                 residual.raw(),
//                                 desired_precision,
//                                 0
//                             );
//                             residual.normalize( mass );
//                             
//                             LOG << "\t\t\t b:       " << rhs_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ab:      " << ( mat * rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t r:       " << residual.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ar:      " << ( mat * residual ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t Ar:      " << ( mat * ( rhs - mat * sol ) ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
//                         }
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
