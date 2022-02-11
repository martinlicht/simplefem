

/**/

#include <iostream>
#include <fstream>
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
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examplesUK.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main()
{
        
        LOG << "Unit Test for Solution of Neumann Problem" << endl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Case 2D" << endl;
            
            LOG << "Initial mesh..." << endl;
            
            MeshSimplicial2D M = UnitedKingdom();

            M.check();
            
            LOG << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_rhs;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_grad;
            std::vector<std::function<FloatVector(const FloatVector&)>> experiments_sol;


            
            // std::function<FloatVector(const FloatVector&) scalarfield = 
            
            Float xfeq = 1.;
            Float yfeq = 1.;
            

            experiments_sol.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) });
                }
            );

            experiments_grad.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ),
                            -yfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] ), 
                        });
                }
            );

            experiments_rhs.push_back( 
                [xfeq,yfeq](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                        +
                        yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] )
                     });
                }
            );

            

            assert( experiments_sol.size() == experiments_rhs.size() );

            LOG << "Solving Poisson Problem with Neumann boundary conditions" << endl;

            int max_l = 5;
            int r = 1;

            for( int l = 0; l <= max_l; l++ ){
                
                LOG << "Level: " << l << std::endl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                
                LOG << "...assemble scalar mass matrices" << endl;
        
                SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix scalar_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 0, r );
                
                LOG << "...assemble vector mass matrix" << endl;
        
                SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                
                SparseMatrix vector_massmatrix_fac = FEECBrokenMassMatrixRightFactor( M, M.getinnerdimension(), 1, r-1 );
                
                LOG << "...assemble differential matrix and transpose" << endl;

                SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                LOG << "...assemble inclusion matrix and transpose" << endl;
        
                SparseMatrix incmatrix = LagrangeInclusionMatrix( M, M.getinnerdimension(), r );

                SparseMatrix incmatrix_t = incmatrix.getTranspose();

                LOG << "...assemble stiffness matrix" << endl;
        
                // ProductOperator 
                // auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                // auto op1 = incmatrix_t * diffmatrix_t;
                // auto op2 = op1 * vector_massmatrix;
                // auto op3 = op2 * diffmatrix;
                // auto stiffness = op3 * incmatrix;

                auto opr1 = diffmatrix & incmatrix;
                auto opr  = vector_massmatrix_fac & opr1;
                auto opl  = opr.getTranspose(); 
                auto stiffness = opl & opr;
                
                stiffness.sortentries();
                auto stiffness_csr = MatrixCSR( stiffness );
                
                //auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
                auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                LOG << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << std::endl;

                const auto& function_sol = experiments_sol[0];
                const auto& function_grad= experiments_grad[0];
                const auto& function_rhs = experiments_rhs[0];
                
                LOG << "...interpolate explicit solution and rhs" << endl;
    
                FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                
                FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                
                LOG << "...measure kernel component: " << std::flush;
    
                Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                
                LOG << average_sol << space << average_rhs << endl;

                LOG << "...measure interpolation commutativity" << endl;
    
                Float commutatorerror = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * interpol_sol ) ).norm();
                LOG << "commutator error: " << commutatorerror << endl;
                
                LOG << "...compute norms of solution and right-hand side:" << endl;
    
                Float sol_norm = ( scalar_massmatrix_fac * interpol_sol ).norm();
                Float rhs_norm = ( scalar_massmatrix_fac * interpol_rhs ).norm();
                
                LOG << "solution norm: " << sol_norm << endl;
                LOG << "rhs norm:      " << rhs_norm << endl;

                LOG << "...create RHS vector" << endl;

                FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                FloatVector sol( M.count_simplices(0), 0. );
                
                LOG << "...iterative solver" << endl;
                
                if(false){
                    sol.zero();
                    timestamp start = gettimestamp();
                    ConjugateResidualMethod CRM( stiffness_csr );
                    CRM.print_modulo = 1+sol.getdimension()/1000;
                    CRM.threshold = 1e-10;
                    CRM.solve( sol, rhs );
                    timestamp end = gettimestamp();
                    LOG << "\t\t\t " << timestamp2measurement( end - start ) << std::endl;
                }
                        
                if(false)
                {
                    sol.zero();
                    timestamp start = gettimestamp();
                    PreconditionedConjugateResidualMethod PCRM( stiffness_csr, stiffness_invprecon );
                    PCRM.print_modulo = 1+sol.getdimension()/10;
                    PCRM.threshold = 1e-10;
                    PCRM.solve( sol, rhs );
                    timestamp end = gettimestamp();
                    LOG << "\t\t\t " << timestamp2measurement( end - start ) << std::endl;
                }

                LOG << "...compute error and residual:" << endl;

                Float errornorm     = ( scalar_massmatrix_fac * ( interpol_sol  - incmatrix * sol ) ).norm();
                Float graderrornorm = ( vector_massmatrix_fac * ( interpol_grad - diffmatrix * incmatrix * sol ) ).norm();
                Float residualnorm  = ( rhs - stiffness * sol ).norm();
                
                // FloatVector gradfoo = diffmatrix * ( interpol_sol - incmatrix * sol );
                // Float graderrornorm = gradfoo.scalarproductwith( vector_massmatrix * gradfoo );
                // Float errornorm1 = interpol_sol * ( scalar_massmatrix * interpol_sol );
                // Float errornorm2 = power_numerical( ( scalar_massmatrix_fac * interpol_sol ).norm(), 2. );

                LOG << "error:     " << errornorm    << endl;
                LOG << "graderror: " << graderrornorm << endl;
                LOG << "residual:  " << residualnorm << endl;


                {
            
                    fstream fs( adaptfilename("./afempoissonneumannUK.vtk"), std::fstream::out );
        
                    VTKWriter vtk( M, fs, "Poisson-Neumann problem" );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();

                    // vtk.write
                    
                    vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                    
                    fs.close();
            
                }
                
                LOG << "Refinement..." << endl;
            
                if( l != max_l ) {

                    FloatVector vec = interpol_grad;

                    FloatVector cellwisemass = 
                        FEECBrokenMassMatrix_cellwisemass( M, M.getinnerdimension(), 1, 0, vec )
                        +
                        Interpolation( M, M.getinnerdimension(), 0, 0, function_rhs );

                    Float maxcellwisemass = cellwisemass.maxnorm();

                    std::vector<int> marked_edges;
                    marked_edges.reserve( 3 * M.count_edges() );

                    for( int s = 0; s < M.count_triangles(); s++ ) 
                    if( cellwisemass.at(s) > 0.5 * maxcellwisemass )
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
        
        
        
        
        LOG << "Finished Unit Test" << endl;
        
        return 0;
}
