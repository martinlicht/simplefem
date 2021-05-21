

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
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        cout << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << endl;
        
        cout << std::setprecision(10);

        if(true){

            cout << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            cout << "Prepare scalar fields for testing..." << endl;
            
            
            const Float xfeq = 1.;
            const Float yfeq = 1.;

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        
                        xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                        +
                        yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                     });
                };
            

            
            
            

            cout << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            ConvergenceTable contable;
            

            int min_l = 2; int max_l = 9;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    cout << "...assemble matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    scalar_massmatrix.sortandcompressentries();
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    vector_massmatrix.sortandcompressentries();
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    auto composed_stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    auto composed_mass      = incmatrix_t * scalar_massmatrix * incmatrix;


                    /*
                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness_prelim = opl & ( vector_massmatrix & opr );
                    stiffness_prelim.sortentries();
                    */
                    //auto stiffness = MatrixCSR( stiffness_prelim );
                    const auto& stiffness = composed_stiffness; 
                    
                    {

                        FloatVector sol_original( M.count_simplices(0), 0. );
                        sol_original.random();
                        sol_original.normalize( composed_mass );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs_original = incmatrix_t * ( scalar_massmatrix * interpol_rhs );
                        rhs_original.zero();
                        
                        if(false)
                        {
                            cout << "CGM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            FloatVector rhs = rhs_original;
                            ConjugateGradientMethod Solver( stiffness );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << sol.norm( composed_mass );
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );
                        }

                        if(false)
                        {
                            cout << "CRM C++" << endl;
                        
                            FloatVector sol = sol_original;
                            FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.print_modulo        = 1;
                            Solver.threshold        = 10000 * machine_epsilon;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve_robust( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << sol.norm( composed_mass );
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );
                        }

                        {
                            cout << "MINRES C++" << endl;
                        
                            FloatVector sol = sol_original;
                            FloatVector rhs = rhs_original;
                            MinimumResidualMethod Solver( stiffness );
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.print_modulo        = 1;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << sol.norm( composed_mass );
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );
                        }

                        if(false)
                        {
                            cout << "HERZOG SOODHALTER C++" << endl;
                        
                            FloatVector sol = sol_original;
                            FloatVector rhs = rhs_original;
                            HerzogSoodhalterMethod Solver( stiffness );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << sol.norm( composed_mass );
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );
                        }
                        
                        
                        contable << nl;
                        
                        contable.print( std::cout, false );

                    }
                    
                }

                cout << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();

                
            } 
        
        }
        
        
        
        
        cout << "Finished Unit Test" << endl;
        
        return 0;
}
