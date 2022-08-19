

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
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
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D Mx = StandardSquare2D();
            
            Mx.check();
            
//             Mx.automatic_dirichlet_flags();

            
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

            
            ConvergenceTable contable("Mass error");
            
            contable << "#nullvec";
                        

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
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << "/" << max_r << std::endl;
                    
                    LOG << "...assemble matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    LOG << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    LOG << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble stiffness matrix" << endl;
            
                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    
//                     auto stiffness_prelim = opl & ( vector_massmatrix & opr );
//                     stiffness_prelim.sortentries();
                    auto stiffness = MatrixCSR( opl & ( vector_massmatrix & opr ) ); // MatrixCSR( stiffness_prelim );
                    
                    auto physical_mass = MatrixCSR( incmatrix_t & ( scalar_massmatrix & incmatrix ) );
                    
                    auto stiffness_diagonal = SparseMatrix( DiagonalOperator( vector_massmatrix.diagonal() ) );
                    assert( stiffness_diagonal.issquare() );
                    assert( stiffness_diagonal.getdimin() == opr.getdimout() );
                    auto simplified_stiffness = MatrixCSR( opl & ( stiffness_diagonal & opr ) );
                    
                    auto idea_prelim = opl & opr;
                    idea_prelim.sortentries();
                    auto idea = MatrixCSR( idea_prelim );
                    
                    const auto& SystemMatrix = stiffness;
//                     const auto& SystemMatrix = simplified_stiffness;
                    
                    
                    assert( SystemMatrix.diagonal().isfinite() );
                    
//                     LOG << SystemMatrix << endl;
                    
                    std::vector<FloatVector> nullvectorgallery;
                    
                    const auto& mass = physical_mass;
//                     const auto& mass = IdentityMatrix(physical_mass.getdimin());
                    
                    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
                    {
                        
                        FloatVector candidate( opr.getdimin(), 0. ); 
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
                            FloatVector rhs( opr.getdimin(), 0. );
                        
                            FloatVector residual( rhs );
                            
                            for( int t = 0; t < max_number_of_purifications; t++ )
                            {
                                
                                ConjugateResidualSolverCSR( 
                                    candidate.getdimension(), 
                                    candidate.raw(), 
                                    rhs.raw(), 
                                    SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
                                    residual.raw(),
                                    desired_precision,
                                    0
                                );
                                
                                LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm() << std::endl;
                                LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + SystemMatrix * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) res:       " << residual.norm() << std::endl;
                                LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) x:         " << candidate.norm() << std::endl;
                                LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) Ax:        " << ( SystemMatrix * candidate ).norm() << std::endl;
                                LOG << "\t\t\t (mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << std::endl;
                                LOG << "\t\t\t (eucl) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm() << std::endl;
                                LOG << "\t\t\t (mass) b - Ax:    " << ( SystemMatrix * candidate - rhs ).norm( mass ) << std::endl;
                                
                                
                                candidate.normalize( mass );

                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << std::endl;
                                LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << std::endl;
                                LOG << "\t\t\t (norm eucl) Ax:        " << ( SystemMatrix* candidate ).norm() << std::endl;
                                LOG << "\t\t\t (norm mass) Ax:        " << ( SystemMatrix * candidate ).norm( mass ) << std::endl;
                                
                                
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
//                                 candidate.normalize( mass );
                                
                                
                                
                                
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
                        
                        if( residual_mass > 1e-6 ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << std::endl;
                            continue;
                        }
                        
                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << std::endl;
                    
                        
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
                        
                        vtk.writeVertexScalarData( nullvector,  "nullvector" , 1.0 );
                        
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
//                                 SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                 residual.raw(),
//                                 desired_precision,
//                                 1
//                             );
//                             sol.normalize( mass );
//                             
//                             assert( sol.isfinite() );
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( SystemMatrix * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( SystemMatrix * sol - rhs ).norm( mass ) << std::endl;
//                         
//                         }
//                         
//                         
//                         
//                         contable << sol.norm( mass ) << ( SystemMatrix * sol ).norm( mass );
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
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
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
//                                 SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                 residual.raw(),
//                                 machine_epsilon,
//                                 1
//                             );
//                             sol.normalize( mass );
//                             
//                             assert( sol.isfinite() );
//                             
//                             LOG << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax_0:      " << ( SystemMatrix * sol_original ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax_0:  " << ( SystemMatrix * sol_original - rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( SystemMatrix * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( SystemMatrix * sol - rhs ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( SystemMatrix * sol ).norm( mass );
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
//                                 SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
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
//                                 SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                 residual.raw(),
//                                 machine_epsilon,
//                                 0
//                             );
//                             sol.normalize( mass );
//                             
//                             LOG << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax_0:      " << ( SystemMatrix * sol_original ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax_0:  " << ( SystemMatrix * sol_original - rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ax:        " << ( SystemMatrix * sol ).norm( mass ) << std::endl;
//                             LOG << "\t\t\t b - Ax:    " << ( SystemMatrix * sol - rhs ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( SystemMatrix * sol ).norm( mass );
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
//                                 SystemMatrix.getA(), SystemMatrix.getC(), SystemMatrix.getV(),
//                                 residual.raw(),
//                                 desired_precision,
//                                 0
//                             );
//                             residual.normalize( mass );
//                             
//                             LOG << "\t\t\t b:       " << rhs_original.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ab:      " << ( SystemMatrix * rhs ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t r:       " << residual.norm( mass ) << std::endl;
//                             LOG << "\t\t\t Ar:      " << ( SystemMatrix * residual ).norm( mass ) << std::endl;
//                             
//                             LOG << "\t\t\t Ar:      " << ( SystemMatrix * ( rhs - SystemMatrix * sol ) ).norm( mass ) << std::endl;
//                             
//                             contable << sol.norm( mass ) << ( SystemMatrix * sol ).norm( mass );
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
