

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
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test for Solution of Dirichlet Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = LShapedDomain2D();

            M.automatic_dirichlet_flags();
            
            M.check();
            
            M.check_dirichlet_flags();
            
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            // The solution of Laplacian problems over L-shaped domainswith a singular function boundary integral method
            // https://onlinelibrary.wiley.com/doi/pdf/10.1002/cnm.489?casa_token=KTbdSboKSK8AAAAA:ISbMXTrwR6i-CocYB6hgQdxdbGgjQxo1QMxRA-L97XFrW_BuEiyUxnXZSVM_SF3DTLmHyGe0ZNdLZtR3
                    
            Float xfeq = 1.;
            Float yfeq = 1.;
            


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
            

            

            

            cout << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            int min_l = 2; 
            int max_l = 5;
            
            ConvergenceTable contable;
            
            const int r = 1;
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                cout << "...assemble matrices" << endl;
        
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
                
                cout << "...iterative solver 1" << endl;
                
                {
                    sol.zero();
                    MinimumResidualMethod Solver( stiffness_csr );
                    Solver.print_modulo        = 1+sol.getdimension();
                    Solver.max_iteration_count = 4 * sol.getdimension();
                    timestamp start = gettimestamp();
                    Solver.solve( sol, rhs );
                    timestamp end = gettimestamp();
                    std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                }

                cout << "...iterative solver 2" << endl;
                
                {
                    aug_sol.zero();
                    MinimumResidualMethod Solver( aug_stiffness_csr );
                    Solver.print_modulo        = 1+aug_sol.getdimension();
                    Solver.max_iteration_count = 4 * aug_sol.getdimension();
                    timestamp start = gettimestamp();
                    Solver.solve( aug_sol, aug_rhs );
                    timestamp end = gettimestamp();
                    std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                }

                cout << "...compute error and residual:" << endl;

                FloatVector error     = aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol;
                FloatVector graderror = aug_diffmatrix * ( aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol );
                Float errornorm       = std::sqrt( error * ( aug_scalar_massmatrix * error ) );
                Float graderrornorm   = std::sqrt( graderror * ( aug_vector_massmatrix * graderror ) );
                
                cout << "error:     " << errornorm     << endl;
                cout << "graderror: " << graderrornorm << endl;
                
                
                        
                        
                contable << errornorm << graderrornorm << nl;
                
                contable.print( std::cout );


                if( r == 1 ){
            
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
        
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    
                    vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                    
                    fs.close();
            
                }
                
                cout << "Refinement..." << endl;
                
                M.uniformrefinement();
                

            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
