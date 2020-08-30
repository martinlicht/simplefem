

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
#include "../../vtk/vtkwriter.mesh2D.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/cgm.hpp"
#include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
#include "../../solver/minres.hpp"
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
            
            MeshSimplicial2D M = UnitSquare2D();
            
            M.check();
                        
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            Float xfeq = 1.;
            Float yfeq = 1.;
            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ) });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector( { 
                            xfeq * Constants::twopi * std::cos( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] ),
                            yfeq * Constants::twopi * std::sin( xfeq * Constants::twopi * vec[0] ) * std::cos( yfeq * Constants::twopi * vec[1] ), 
                        });
                };
            

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
            

            int min_l = 4; int max_l = 4;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    cout << "...assemble matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    cout << "...assemble vector mass matrix" << endl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    cout << "...assemble differential matrix and transpose" << endl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    cout << "...assemble inclusion matrix and transpose" << endl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble stiffness matrix" << endl;
            
                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    
                    auto stiffness_prelim = opl & ( vector_massmatrix & opr );
                    stiffness_prelim.sortentries();
                    auto stiffness = MatrixCSR( stiffness_prelim );
                    
                    auto idea_prelim = opl & opr;
                    idea_prelim.sortentries();
                    auto idea = MatrixCSR( idea_prelim );
                    
                    const auto& mat = stiffness;
                    
                    {

                        FloatVector sol_original( opr.getdimin(), 0. );
                        FloatVector rhs_original( opr.getdimin(), 0. );
                        
                        sol_original.random(); sol_original.normalize();
                        rhs_original.random(); rhs_original.normalize();
                        
                        {
                            cout << "Filter out from x (CGM)" << endl;
                        
                            FloatVector sol( sol_original );
                            FloatVector rhs( rhs_original.getdimension(), 0. );
                            FloatVector residual( rhs );
                            
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mat.getA(), mat.getC(), mat.getV(),
                                residual.raw(),
                                machine_epsilon,
                                0
                            );
                            sol.normalize();
                            
                            std::cout << "\t\t\t x_0:       " << sol_original.norm() << std::endl;
                            std::cout << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm() << std::endl;
                            std::cout << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm() << std::endl;
                            
                            std::cout << "\t\t\t x:         " << sol.norm() << std::endl;
                            std::cout << "\t\t\t Ax:        " << ( mat * sol ).norm() << std::endl;
                            std::cout << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm() << std::endl;
                            
                            contable << sol.norm() << ( mat * sol ).norm();
                        }

                        {
                            cout << "Filter out from x (CRM)" << endl;
                        
                            FloatVector sol( sol_original );
                            FloatVector rhs( rhs_original.getdimension(), 0. );
                            FloatVector residual( rhs );
                            
                            ConjugateResidualSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mat.getA(), mat.getC(), mat.getV(),
                                residual.raw(),
                                machine_epsilon,
                                0
                            );
                            sol.normalize();
                            
                            std::cout << "\t\t\t x_0:       " << sol_original.norm() << std::endl;
                            std::cout << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm() << std::endl;
                            std::cout << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm() << std::endl;
                            
                            std::cout << "\t\t\t x:         " << sol.norm() << std::endl;
                            std::cout << "\t\t\t Ax:        " << ( mat * sol ).norm() << std::endl;
                            std::cout << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm() << std::endl;
                            
                            contable << sol.norm() << ( mat * sol ).norm();
                        }

                        {
                            cout << "Filter out from b" << endl;
                        
                            FloatVector sol( sol_original.getdimension(), 0. );
                            FloatVector rhs( rhs_original );
                            FloatVector residual( rhs );
                            
                            ConjugateResidualSolverCSR_textbook( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mat.getA(), mat.getC(), mat.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            residual.normalize();
                            
                            std::cout << "\t\t\t b:       " << rhs_original.norm() << std::endl;
                            std::cout << "\t\t\t Ab:      " << ( mat * rhs ).norm() << std::endl;
                            
                            std::cout << "\t\t\t r:       " << residual.norm() << std::endl;
                            std::cout << "\t\t\t Ar:      " << ( mat * residual ).norm() << std::endl;
                            
                            std::cout << "\t\t\t Ar:      " << ( mat * ( rhs - mat * sol ) ).norm() << std::endl;
                            
                            contable << sol.norm() << ( mat * sol ).norm();
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
