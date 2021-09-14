

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
#include "../../solver/chebyshev.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Compare numerical solvers CRM vs MINRES for Solution of Dirichlet Problem" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
        
        LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            LOG << "Prepare scalar fields for testing..." << endl;
            

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
            

            
            
            

            LOG << "Solving Poisson Problem with Dirichlet boundary conditions" << endl;

            ConvergenceTable contable;
            

            int min_l = 7; int max_l = 7;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
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
                    auto stiffness_prelim = opl & ( vector_massmatrix & opr );
                    stiffness_prelim.sortentries();
                    auto stiffness = MatrixCSR( stiffness_prelim );
                    
                    {

                        FloatVector sol( M.count_simplices(0), 0. );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        if(false)
                        {
                            LOG << "CGM - CSR Classic" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        if(false)
                        {
                            LOG << "CGM - CSR Classic with diagonal preconditioning" << endl;
                            
                            auto precon = InverseDiagonalPreconditioner( stiffness );

                            
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                precon.getdiagonal().raw()
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        {
                            LOG << "CGM - CSR Classic with SSOR" << endl;
                            
                            auto diagonal = stiffness.diagonal();

                            
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR_SSOR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                1.0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        if(false)
                        {
                            LOG << "CRM - CSR Classic" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateResidualSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        if(false)
                        {
                            LOG << "CRM - CSR Textbook" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateResidualSolverCSR_textbook( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        if(false)
                        {
                            LOG << "MINRES CSR" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            MINRESCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }


                        if(false)
                        {
                            LOG << "WHATEVER CSR" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            WHATEVER( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                stiffness.getA(), stiffness.getC(), stiffness.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( stiffness * sol - rhs ).norm() );;
                        }

                        
                        
                        contable << nl;
                        
                        contable.lg( false );

                    }
                    
                }

                LOG << "Refinement..." << endl;
            
                if( l != max_l ) M.uniformrefinement();

                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
