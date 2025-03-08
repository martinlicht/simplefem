

/**/

#include <cmath>

#include <algorithm>
#include <functional>
#include <vector>

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/densematrix.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test for Solution of Neumann Problem" << nl;
    
    if(true){

        LOG << "Case 2D" << nl;
        
        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();

        for( int t = 0; t < 3; t++ ) M.uniformrefinement();
        
        M.check();
        
        LOG << "Prepare scalar fields for testing..." << nl;
        

        std::function<FloatVector(const FloatVector&)> constant_one
            = [](const FloatVector& vec) -> FloatVector {
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 1. });
                };
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_rhs;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_sol;


        
        Float xfeq = 1.;
        Float yfeq = 1.;
        

        experiments_sol.push_back( 
            [xfeq,yfeq](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 2 );
                // return FloatVector({ 1. });
                return FloatVector({ std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) });
            }
        );

        experiments_grad.push_back( 
            [xfeq,yfeq](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 2 );
                return FloatVector( { 
                        -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ),
                        -yfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] ), 
                    });
            }
        );

        experiments_rhs.push_back( 
            [xfeq,yfeq](const FloatVector& vec) -> FloatVector {
                assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                    +
                    yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                    });
            }
        );

        

        assert( experiments_sol.size() == experiments_rhs.size() );

        const int min_l = 1; 

        const int max_l = 5;
        
        const int r = 2;

        for( int l = 0; l < min_l; l++ )
            M.uniformrefinement();

        for( int l = min_l; l <= max_l; l++ ){
            
            LOG << "Level: " << l << "/" << max_l << nl;
            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
            
            
            LOG << "...assemble scalar mass matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
            
            SparseMatrix scalar_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
            
            LOG << "...assemble vector mass matrix" << nl;
    
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
            
            SparseMatrix vector_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 1, r-1 );
            
            LOG << "...assemble differential matrix and transpose" << nl;

            SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

            SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

            LOG << "...assemble inclusion matrix and transpose" << nl;
    
            SparseMatrix incmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );

            SparseMatrix incmatrix_t = incmatrix.getTranspose();

            LOG << "...assemble stiffness matrix" << nl;
    
            auto opr1 = diffmatrix & incmatrix;
            auto opr  = vector_massmatrix_fac & opr1;
            auto opl  = opr.getTranspose(); 
            auto stiffness = opl & opr;
            
            stiffness.sortentries();
            auto stiffness_csr = MatrixCSR( stiffness );
            
            //auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
            auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
            LOG << "Average value of diagonal preconditioner: " << stiffness_invprecon.getDiagonal().average() << nl;

            const auto& function_sol = experiments_sol[0];
            const auto& function_grad= experiments_grad[0];
            const auto& function_rhs = experiments_rhs[0];
            
            LOG << "...interpolate explicit solution and rhs" << nl;

            FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
            FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
            FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
            
            FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
            
            LOG << "...measure kernel component: ";

            Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
            Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
            
            LOG << average_sol << space << average_rhs << nl;

            LOG << "...create RHS vector" << nl;

            FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

            FloatVector sol( M.count_simplices(0), 0. );
            
            LOG << "...iterative solver" << nl;
            
            {
                sol.zero();
                timestamp start = timestampnow();
                PreconditionedConjugateResidualMethod PCRM( stiffness_csr, stiffness_invprecon );
                PCRM.print_modulo = 1+sol.getdimension()/10;
                PCRM.tolerance = desired_precision;
                PCRM.solve( sol, rhs );
                timestamp end = timestampnow();
                LOG << "\t\t\t " << timestamp2measurement( end - start ) << nl;
            }

            {
                sol.zero();
                timestamp start = timestampnow();
                ConjugateResidualMethod CRM( stiffness_csr );
                CRM.print_modulo = 1+sol.getdimension()/1000;
                CRM.tolerance = desired_precision;
                CRM.solve( sol, rhs );
                timestamp end = timestampnow();
                LOG << "\t\t\t " << timestamp2measurement( end - start ) << nl;
            }
                    
            LOG << "...compute error and residual" << nl;

            Float errornorm     = ( scalar_massmatrix_fac * ( interpol_sol  - incmatrix * sol ) ).norm();
            Float graderrornorm = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * incmatrix * sol ) ).norm();
            Float residualnorm  = ( rhs - stiffness * sol ).norm();
            
            LOG << "error:     " << errornorm    << nl;
            LOG << "graderror: " << graderrornorm << nl;
            LOG << "residual:  " << residualnorm << nl;


            if( l != max_l ) {

                LOG << "Refinement..." << nl;
        
                FloatVector vec = interpol_grad - diffmatrix * ( incmatrix * sol );

                FloatVector cellwisemass = 
                    FEECBrokenMassMatrix_cellwisemass( M, M.getinnerdimension(), 1, 0, vec )
                    +
                    Interpolation( M, M.getinnerdimension(), 0, 0, function_rhs );

                Float maxcellwisemass = cellwisemass.maxnorm();

                std::vector<int> marked_edges;
                marked_edges.reserve( 3 * M.count_edges() );

                for( int s = 0; s < M.count_triangles(); s++ ) 
                if( cellwisemass.at(s) > 0.1 * maxcellwisemass )
                {
                    // LOG << M.get_triangle_edge( s, 0 ) << space << M.get_triangle_edge( s, 1 ) << space << M.get_triangle_edge( s, 2 ) << nl;
                    marked_edges.push_back( M.get_triangle_edge( s, 0 ) );
                    marked_edges.push_back( M.get_triangle_edge( s, 1 ) );
                    marked_edges.push_back( M.get_triangle_edge( s, 2 ) );
                }

                std::sort( marked_edges.begin(), marked_edges.end() );
                auto temp = std::unique( marked_edges.begin(), marked_edges.end() );
                marked_edges.erase( temp, marked_edges.end() );
                
                LOG << "marked edges: " << marked_edges.size() << "/" << M.count_edges() << nl;

                // for( int e = 0; e < M.count_edges(); e++ )
                //     if( e % 10 == 0)
                //         marked_edges.push_back( e );
                
                M.newest_vertex_bisection_recursive( marked_edges );

            }
            
        } 
    
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
