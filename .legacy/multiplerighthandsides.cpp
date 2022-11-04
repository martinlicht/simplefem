

/**/

#include <cmath>
// #include <ostream>
// #include <fstream>
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
// #include "../../solver/chebyshev.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            ConvergenceTable contable;
            
            contable << "Time" << "L2 Error" << "L2 Residual" << "Iteration percentage" << "Angle" << "Ratio";

            contable.display_convergence_rates = false;
            
            const int min_l = 0;
            
            const int max_l = 0;

            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    LOG << "...assemble matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    LOG << "...assemble differential matrix and transpose" << nl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << nl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble stiffness and mass matrices" << nl;
            
                    const auto composed_stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    const auto composed_mass      = incmatrix_t * scalar_massmatrix * incmatrix;

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness_csr_prelim = opl & ( vector_massmatrix & opr );
                    stiffness_csr_prelim.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness_csr_prelim );

                    const auto& stiffness = stiffness_csr;
                    const auto& mass      = composed_mass;
                    
                    {
                        
                        int number_of_problems = 20;

                        vector<FloatVector> sol_original( number_of_problems, FloatVector( M.count_simplices(0), 0. ) );
                        vector<FloatVector> rhs_original( number_of_problems, FloatVector( M.count_simplices(0), 0. ) );

                        for( int i = 0; i < number_of_problems; i++ ) {
                            sol_original[i].random();
                            rhs_original[i] = stiffness * sol_original[i];
                            Float norm = rhs_original[i].norm();
                            rhs_original[i] /= norm;
                            sol_original[i] /= norm;
                        }

                        DenseMatrix Angles( number_of_problems );
                        for( int i = 0; i < number_of_problems; i++ )
                        for( int j = 0; j < number_of_problems; j++ )
                        {
                            Float v = acos( rhs_original[i] * rhs_original[j] ) / Constants::pi;
                            Float w = ( rhs_original[i] - rhs_original[j] ).norm();
                            Angles(i,j) = v;
                        }
                        LOG << Angles.text();

// Use Prim's algorithm to grow a minimum weight spanning tree 
// Start with one vertex; given a tree, add the minumum weight edge to that connects the tree to a new vertex
// can be implemented with some matrix 

                        
                        {
                            LOG << "CGM - CSR Classic with SSOR" << nl;
                            
                            const auto diagonal = stiffness.diagonal();
                            FloatVector residual( M.count_simplices(0), 0. );
                            auto max_iteration_count = M.count_simplices(0);
                            
                            for( int i = 0; i < number_of_problems; i++ )
                            {
                                
                                const FloatVector prev_rhs = (i!=0) ? rhs_original[i-1] : FloatVector( M.count_simplices(0), 0. );
                                const FloatVector curr_rhs = rhs_original[i];

                                const Float rhs_angle = acos( curr_rhs * prev_rhs ) / Constants::pi;
                                
                                const FloatVector prev_sol = (i!=0) ? sol_original[i-1] : FloatVector( M.count_simplices(0), 0. );
                                FloatVector       curr_sol = prev_sol;// FloatVector( M.count_simplices(0), 0. );

                                
                                timestamp start = gettimestamp();
                                auto recent_iteration_count = 
                                ConjugateGradientSolverCSR_SSOR( 
                                    curr_sol.getdimension(), 
                                    curr_sol.raw(), 
                                    curr_rhs.raw(), 
                                    stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                    residual.raw(),
                                    desired_precision,
                                    -1,
                                    diagonal.raw(),
                                    1.0
                                );

                                timestamp end = gettimestamp();
                                // LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                                
                                auto stat_tmp  = static_cast<Float>( end - start );
                                
                                auto stat_sol = Float( ( curr_sol - sol_original[i] ).norm() );
                                
                                auto stat_res = Float( ( stiffness * curr_sol - curr_rhs ).norm() );
                                
                                auto stat_num = Float( recent_iteration_count ) / max_iteration_count;
                                
                                auto stat_rat = rhs_angle / recent_iteration_count;
                            
                                contable << stat_tmp << stat_sol << stat_res << stat_num << rhs_angle << stat_rat << nl;

                            }

                        contable.lg();

                        }
                    
                    }
                    
                }

                if( l != max_l ){ 
                    LOG << "Refinement..." << nl;
                    M.uniformrefinement();
                }

                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test" << nl;
        
        return 0;
}
