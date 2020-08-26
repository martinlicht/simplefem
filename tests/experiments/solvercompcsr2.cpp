

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
#include "../../solver/chebyshev.hpp"
#include "../../solver/cgm.hpp"
#include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
#include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.lagrangeincl.hpp"
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
            

            int min_l = 2; int max_l = 10;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                cout << "Level: " << l << std::endl;
                cout << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    cout << "...assemble matrices" << endl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix incmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    cout << "...assemble global mass matrix" << endl;
            

                    auto mass_prelim = incmatrix_t & ( scalar_massmatrix & incmatrix );
                    mass_prelim.sortentries();
                    auto mass = MatrixCSR( mass_prelim );
                    
                    {

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        {
                            cout << "CGM - CSR Classic" << endl;
                        
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }

                        {
                            cout << "CRM - CSR Classic" << endl;
                        
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            ConjugateResidualSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }

                        {
                            cout << "CRM - CSR Textbook" << endl;
                        
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            ConjugateResidualSolverCSR_textbook( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }

                        {
                            cout << "MINRES CSR" << endl;
                        
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            MINRESCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }


                        {
                            cout << "WHATEVER CSR" << endl;
                        
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            WHATEVER( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }


                        {
                            cout << "CGM diagonal preconditioner CSR" << endl;
                        
                            FloatVector invprecon = mass_prelim.InverseDiagonalPreconditioner();
                            invprecon.setentries( 1. );
                            assert( invprecon.isfinite() );
                            assert( invprecon.ispositive() );
                            
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                0,
                                invprecon.raw()
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start) << Float( ( mass * sol - rhs ).norm() );;
                        }
                        
                        
                        if(false){
                            cout << "CHEBYSHEV CSR" << endl;
                        
                            FloatVector invprecon = mass_prelim.InverseDiagonalPreconditioner();
                            assert( invprecon.isfinite() );
                            assert( invprecon.ispositive() );
                            
                            sol.zero();
                            timestamp start = gettimestamp();
                            FloatVector residual( rhs );
                            CheybyshevIteration_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                1e-16,
                                1,
                                invprecon.raw(),
                                0.,
                                200 * invprecon.maxnorm()
                            );

                            timestamp end = gettimestamp();
                            std::cout << "\t\t\t Time: " << timestamp2string( end - start ) << std::endl;
                            contable << Float(end - start);
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
