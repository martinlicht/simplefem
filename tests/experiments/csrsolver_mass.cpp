

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
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
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
            

            int min_l = 2; int max_l = 9;

            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << std::endl;
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

//                         if(false) 
                        {
                            LOG << "CGM - CSR Classic" << endl;
                        
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
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
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
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
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
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
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
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
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
                        }


//                         if(false)
                        {
                            LOG << "CGM diagonal preconditioner CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( mass_prelim );
//                             invprecon.setentries( 1. );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().isnonnegative() );
                            
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                invprecon.getdiagonal().raw()
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
                        }
                        
                        
//                         if(false)
                        {
                            LOG << "CGM SSOR preconditioner CSR" << endl;
                        
                            FloatVector diagonal = mass.diagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            ConjugateGradientSolverCSR_SSOR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                0.9123456789
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
                        }
                        
                        
                        if(false)
                        {
                            LOG << "CHEBYSHEV CSR" << endl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( mass_prelim );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().ispositive() );
                            
                            sol.zero();
                            FloatVector residual( rhs );
                            timestamp start = gettimestamp();
                            CheybyshevIteration_DiagonalPreconditioner( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mass.getA(), mass.getC(), mass.getV(),
                                residual.raw(),
                                desired_precision,
                                10,
                                invprecon.getdiagonal().raw(),
                                0.,
                                100 * invprecon.getdiagonal().maxnorm()
                            );

                            timestamp end = gettimestamp();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << std::endl;
                            contable << static_cast<Float>( end - start ) << Float( ( mass * sol - rhs ).norm() );
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
