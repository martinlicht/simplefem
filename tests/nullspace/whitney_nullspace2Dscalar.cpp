

/**/

#include <ostream>
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
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
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
            
//             Mx.automatic_dirichlet_flags();

            
            MeshSimplicial2D M;
            
            for( int i = 0; i < 3; i++ )
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

            
            ConvergenceTable contable("Nullvectors found");
            
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
            
                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    
//                     auto stiffness_prelim = opl & ( vector_massmatrix & opr );
//                     stiffness_prelim.sortentries();
                    auto stiffness = MatrixCSR( opl & ( vector_massmatrix & opr ) ); // MatrixCSR( stiffness_prelim );
                    
                    auto physical_mass = MatrixCSR( incmatrix_t & ( scalar_massmatrix & incmatrix ) );
                    
                    auto stiffness_diagonal = SparseMatrix( DiagonalOperator( vector_massmatrix.getDiagonal() ) );
                    assert( stiffness_diagonal.issquare() );
                    assert( stiffness_diagonal.getdimin() == opr.getdimout() );
                    auto simplified_stiffness = MatrixCSR( opl & ( stiffness_diagonal & opr ) );
                    
                    auto idea_prelim = opl & opr;
                    idea_prelim.sortentries();
                    auto idea = MatrixCSR( idea_prelim );
                    
                    const auto& SystemMatrix = stiffness;
//                     const auto& SystemMatrix = simplified_stiffness;
                    
                    
                    assert( SystemMatrix.getDiagonal().isfinite() );
                    
//                     LOG << SystemMatrix << nl;
                    
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
                            LOG << "\t\t\t Preprocessed mass: " << reduced_mass << nl;
                            
                            if( reduced_mass < 1e-6 ) {
                                LOG << "**** The candidate already has very small mass" << nl;
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
                        
                        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;
                    
                        
                        nullvectorgallery.push_back( candidate );
                    }
                    
                    
                    
                    LOG << "How much nullspace are our vectors?" << nl;
                    for( const auto& nullvector : nullvectorgallery ) {
                        Float mass_norm = ( SystemMatrix * nullvector ).norm(mass);
                        assert( is_numerically_small( mass_norm ) );
                        LOG << mass_norm << tab;
                    }
                    LOG << nl;
                    
                    LOG << "How orthonormal are our vectors?" << nl;
                    for( const auto& nullvector1 : nullvectorgallery ) {
                        for( const auto& nullvector2 : nullvectorgallery ) {
                            Float mass_norm = mass * nullvector1 * nullvector2;
                            assert( is_numerically_small( mass_norm ) );
                            LOG << mass_norm << tab;
                        }
                        LOG << nl;
                    }
                    
                    
                    contable << static_cast<Float>(nullvectorgallery.size());     

                    
                    const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );

                    for( const auto& nullvector : nullvectorgallery )
                    {
                
                        fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
            
                        VTKWriter vtk( M, fs, getbasename(__FILE__) );
                        // vtk.writeCoordinateBlock();
                        // vtk.writeTopDimensionalCells();
                        
                        auto reduced_nullvector = interpol_matrix * incmatrix * nullvector;

                        vtk.writeVertexScalarData( reduced_nullvector,  "nullvector" , 1.0 );
                        
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



