

/**/

#include <fstream>










#include "../basic.hpp"
#include "../utility/files.hpp"
#include "../utility/convergencetable.hpp"
#include "../operators/composedoperators.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../sparse/matcsr.hpp"
#include "../mesh/mesh.simplicial2D.hpp"
#include "../mesh/examples2D.hpp"
#include "../vtk/vtkwriter.hpp"
#include "../solver/sparsesolver.hpp"
#include "../solver/iterativesolver.hpp"
#include "../solver/inv.hpp"
#include "../solver/systemsparsesolver.hpp"
#include "../solver/systemsolver.hpp"
#include "../fem/global.contraction.hpp"
#include "../fem/global.diffmatrix.hpp"
#include "../fem/global.elevation.hpp"
#include "../fem/global.interpol.hpp"
#include "../fem/global.massmatrix.hpp"
#include "../fem/global.sullivanincl.hpp"
#include "../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: 2D Maxwell System" << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M0 = LShapedDomain2D();
            MeshSimplicial2D M1 = LShapedDomain2D();
            auto& M = M0;

            for( int e = 0; e < M.count_edges(); e++ )
            {
                if( M.count_edge_triangle_parents(e) == 2 ) continue;
                assert( M.count_edge_triangle_parents(e) == 1 );

                int v0 = M.get_edge_vertex( e, 0 );
                int v1 = M.get_edge_vertex( e, 1 );

                auto cs0 = M.getcoordinates().getvectorclone( v0 );
                auto cs1 = M.getcoordinates().getvectorclone( v1 );

                if( cs0[0] == cs1[0] ) {
                    M0.set_flag( 0, v0, SimplexFlagDirichlet );
                    M0.set_flag( 0, v1, SimplexFlagDirichlet );
                    M0.set_flag( 1, e,  SimplexFlagDirichlet );
                    LOG << "x: " << e << space << v0 << space << v1 << nl;
                } else if( cs0[1] == cs1[1] ) {
                    M1.set_flag( 0, v0, SimplexFlagDirichlet );
                    M1.set_flag( 0, v1, SimplexFlagDirichlet );
                    M1.set_flag( 1, e,  SimplexFlagDirichlet );
                    LOG << "y: " << e << space << v0 << space << v1 << nl;
                }
                    
            }

            

            
            std::function<FloatVector(const FloatVector&)> function_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    return FloatVector({
                        -1.0
                        ,
                        +1.0
                     });
                };

            std::function<FloatVector(const FloatVector&)> function_x = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    return FloatVector({
                        -1.0
                        ,
                        +1.0
                     });
                };

            std::function<FloatVector(const FloatVector&)> function_y = 
                [=](const FloatVector& vec) -> FloatVector{
                    Assert( vec.getdimension() == 2 );
                    return FloatVector({
                        -1.0
                        ,
                        +1.0
                     });
                };

                
            

            const int min_l = 0; 
            
            const int max_l = 5;
            
            const int min_r = 1; 
            
            const int max_r = 1;
            

            
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
                    
                    LOG << "... assemble matrices" << nl;
            
                    
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r-1 );
                    
                    SparseMatrix scalar_massmatrix_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix0   = FEECSullivanInclusionMatrix( M0, M0.getinnerdimension(), 0, r );
                    SparseMatrix scalar_incmatrix1   = FEECSullivanInclusionMatrix( M1, M1.getinnerdimension(), 0, r );
                    SparseMatrix scalar_incmatrix0_t = scalar_incmatrix0.getTranspose();
                    SparseMatrix scalar_incmatrix1_t = scalar_incmatrix1.getTranspose();

                    FloatVector dir0 = Interpolation( M, 2, 1, 0, function_x );
                    FloatVector dir1 = Interpolation( M, 2, 1, 0, function_y );

                    SparseMatrix contractionmatrix0   = FEECBrokenContractionMatrix( M, 2, 1, r-1, 1, 0, dir0 );
                    SparseMatrix contractionmatrix1   = FEECBrokenContractionMatrix( M, 2, 1, r-1, 1, 0, dir1 );
                    SparseMatrix contractionmatrix0_t = contractionmatrix0.getTranspose();
                    SparseMatrix contractionmatrix1_t = contractionmatrix1.getTranspose();
                    
                    

                    const auto SystemMatrix1 = 
                        Block2x2Operator( 
                            scalar_incmatrix0_t * scalar_diffmatrix_t * contractionmatrix0_t * scalar_massmatrix * contractionmatrix0 * scalar_diffmatrix * scalar_incmatrix0,
                            scalar_incmatrix0_t * scalar_diffmatrix_t * contractionmatrix0_t * scalar_massmatrix * contractionmatrix1 * scalar_diffmatrix * scalar_incmatrix1,
                            scalar_incmatrix1_t * scalar_diffmatrix_t * contractionmatrix1_t * scalar_massmatrix * contractionmatrix0 * scalar_diffmatrix * scalar_incmatrix0,
                            scalar_incmatrix1_t * scalar_diffmatrix_t * contractionmatrix1_t * scalar_massmatrix * contractionmatrix1 * scalar_diffmatrix * scalar_incmatrix1 
                        );
                    
                    const auto SystemMatrix2 = 
                        Block2x2Operator( 
                             scalar_incmatrix0_t * scalar_diffmatrix_t * contractionmatrix1_t * scalar_massmatrix * contractionmatrix1 * scalar_diffmatrix * scalar_incmatrix0,
                            -scalar_incmatrix0_t * scalar_diffmatrix_t * contractionmatrix1_t * scalar_massmatrix * contractionmatrix0 * scalar_diffmatrix * scalar_incmatrix1,
                             scalar_incmatrix1_t * scalar_diffmatrix_t * contractionmatrix0_t * scalar_massmatrix * contractionmatrix1 * scalar_diffmatrix * scalar_incmatrix0,
                            -scalar_incmatrix1_t * scalar_diffmatrix_t * contractionmatrix0_t * scalar_massmatrix * contractionmatrix0 * scalar_diffmatrix * scalar_incmatrix1
                        );
                    
                    const auto SystemMatrix = SystemMatrix1 + SystemMatrix2;
                    
                    
                    
                    
                    {

                        FloatVector interpol_rhs  = concatFloatVector( -dir0, dir1 );
                        
                        FloatVector rhs = concatFloatVector( 
                                scalar_incmatrix0_t * ( scalar_massmatrix_plus * dir0 ),
                                scalar_incmatrix1_t * ( scalar_massmatrix_plus * dir1 )
                            );

                        FloatVector sol( rhs.getdimension(), 0. );
                        
                        
                        timestamp start = gettimestamp();
                        
                        auto minres = MinimumResidualMethod( SystemMatrix );
                        minres.solve( sol, rhs );

                        timestamp end = gettimestamp();
        
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                        
                        assert( sol.isfinite() );

                        // LOG << "...compute residual:" << nl;
                        // LOG << "residual:   " << residualnorm  << nl;

                        
                        if( r == 1 ){
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            vtk.writeCoordinateBlock();
                            vtk.writeTopDimensionalCells();

                            auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );

                            fs.close();
                        }


                    }
                    
                }

                if( l != max_l ) { 
                    LOG << "Refinement..." << nl; 
                    M0.uniformrefinement(); 
                    M1.uniformrefinement(); 
                }

        
            } 
                    
        }
        
        
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}




