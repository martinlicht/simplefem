

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
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test for Solution of Dirichlet Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = LShapedDomain2D();

            M.automatic_dirichlet_flags();
            
            M.check();
            
            M.check_dirichlet_flags();
            
            LOG << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            // The solution of Laplacian problems over L-shaped domainswith a singular function boundary integral method
            // https://onlinelibrary.wiley.com/doi/pdf/10.1002/cnm.489?casa_token=KTbdSboKSK8AAAAA:ISbMXTrwR6i-CocYB6hgQdxdbGgjQxo1QMxRA-L97XFrW_BuEiyUxnXZSVM_SF3DTLmHyGe0ZNdLZtR3
                    
            const Float xfeq = 1.;
            const Float yfeq = 1.;
            


            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({
                        1.0
//                         xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
//                         +
//                         yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                     });
                };
            

            

            

            LOG << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            const int min_l = 0; 
            const int max_l = 3;
            

            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;

            
            const int r = 1;
            
            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                LOG << "...assemble matrices" << endl;
        
                SparseMatrix     scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+0 );
                SparseMatrix aug_scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
                
                SparseMatrix     vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                SparseMatrix aug_vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-0 );
                
                SparseMatrix     diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+0 );
                SparseMatrix aug_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );

                SparseMatrix     diffmatrix_t =     diffmatrix.getTranspose();
                SparseMatrix aug_diffmatrix_t = aug_diffmatrix.getTranspose();

                SparseMatrix     incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+0 );
                SparseMatrix aug_incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );

                SparseMatrix     incmatrix_t =     incmatrix.getTranspose();
                SparseMatrix aug_incmatrix_t = aug_incmatrix.getTranspose();

                
                SparseMatrix elevation_matrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r+0, 1 );
                

                auto opr  = diffmatrix & incmatrix;
                auto opl  = opr.getTranspose(); 
                auto stiffness = opl & ( vector_massmatrix & opr );                
                stiffness.sortentries();
                auto stiffness_csr = MatrixCSR( stiffness );
                
                auto aug_opr  = aug_diffmatrix & aug_incmatrix;
                auto aug_opl  = aug_opr.getTranspose(); 
                auto aug_stiffness = aug_opl & ( aug_vector_massmatrix & aug_opr );                
                aug_stiffness.sortentries();
                auto aug_stiffness_csr = MatrixCSR( aug_stiffness );
                
                
                const auto& function_rhs  = experiment_rhs;
                FloatVector     interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                FloatVector aug_interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r+1, function_rhs  );

                FloatVector     rhs =     incmatrix_t * (     scalar_massmatrix *     interpol_rhs );
                FloatVector aug_rhs = aug_incmatrix_t * ( aug_scalar_massmatrix * aug_interpol_rhs );

                FloatVector     sol(     incmatrix.getdimin(), 0. );
                FloatVector aug_sol( aug_incmatrix.getdimin(), 0. );
                
                timestamp start = gettimestamp();

                {
                    LOG << "...iterative solver 1" << endl;
                
                    sol.zero();
                    ConjugateGradientMethod Solver( stiffness_csr );
                    Solver.solve( sol, rhs );
                }

                {
                    LOG << "...iterative solver 2" << endl;
                
                    aug_sol.zero();
                    ConjugateGradientMethod Solver( aug_stiffness_csr );
                    Solver.solve( aug_sol, aug_rhs );
                }

                timestamp end = gettimestamp();
                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                LOG << "...compute error and residual:" << endl;

                FloatVector error     = aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol;
                FloatVector graderror = aug_diffmatrix * ( aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol );
                Float errornorm       = std::sqrt( error * ( aug_scalar_massmatrix * error ) );
                Float graderrornorm   = std::sqrt( graderror * ( aug_vector_massmatrix * graderror ) );
                Float residualnorm  = ( rhs - stiffness * sol ).norm();

                LOG << "error:     " << errornorm    << endl;
                LOG << "graderror: " << graderrornorm << endl;
                LOG << "residual:  " << residualnorm << endl;
                LOG << "time:      " << Float( end - start ) << endl;
                        
                contable << errornorm;
                contable << graderrornorm;
                contable << residualnorm;
                contable << Float( end - start );
                contable << nl;
                
                contable.lg();


                if( r == 1 ){
            
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
        
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                    vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                    
                    fs.close();
            
                }
                
                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
                

            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
