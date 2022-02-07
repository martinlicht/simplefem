

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
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/herzogsoodhalter.hpp"
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
        
        LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            
            
            LOG << "Prepare scalar fields for testing..." << endl;
            

            // std::function<FloatVector(const FloatVector&)> constant_one
            //     = [](const FloatVector& vec) -> FloatVector{
            //             assert( vec.getdimension() == 2 );
            //             return FloatVector({ 1. });
            //         };
            
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
            

            
            
            

            LOG << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            // ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");

            // contable_sol.print_transpose_instead_of_standard = true;
            contable_res.print_transpose_instead_of_standard = true;
            contable_num.print_transpose_instead_of_standard = true;
            

            const int min_l = 2;
            
            const int max_l = 7;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    LOG << "...assemble matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble global mass matrix" << endl;
            

                    auto mass_prelim = incmatrix_t & ( scalar_massmatrix & incmatrix );
                    mass_prelim.sortentries();
                    auto mass = MatrixCSR( mass_prelim );
                    
                    {

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        {
                            LOG << "CGM C++" << endl;
                        
                            sol.zero();
                            ConjugateGradientMethod Solver( mass );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << sol.norm( mass );

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( mass * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        {
                            LOG << "CRM C++" << endl;
                        
                            sol.zero();
                            ConjugateResidualMethod Solver( mass );
                            Solver.print_modulo        = 0;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            
                            LOG << sol.norm( mass );

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( mass * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        if(false)
                        {
                            LOG << "MINRES C++" << endl;
                        
                            sol.zero();
                            MinimumResidualMethod Solver( mass );
                            Solver.print_modulo        = 1;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << sol.norm( mass );

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( mass * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }

                        {
                            LOG << "HERZOG SOODHALTER C++" << endl;
                        
                            sol.zero();
                            HerzogSoodhalterMethod Solver( mass );
                            Solver.print_modulo        = 1;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * sol.getdimension();
                            timestamp start = gettimestamp();
                            Solver.solve( sol, rhs );
                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;

                            LOG << sol.norm( mass );

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( mass * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                        }
                        
                        
                        // contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                    
                        // contable_sol.lg( false );
                        contable_res.lg( false );
                        contable_num.lg( false );

                    }
                    
                }

                if( l != max_l ){ 
                    LOG << "Refinement..." << endl;
                    M.uniformrefinement();
                }

                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
