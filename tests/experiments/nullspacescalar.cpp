

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
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
//             {
//                 auto M2 = M;
//                 M2.getcoordinates().shift( { 3.0, 0.0 } );
//                 M.merge( M2 );
//             }
                        
            cout << "Prepare scalar fields for testing..." << endl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            

            cout << "Nullspace computation" << endl;

            ConvergenceTable contable;
            

            int min_l = 0; int max_l = 0;

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
                    
//                     auto stiffness_prelim = opl & ( vector_massmatrix & opr );
//                     stiffness_prelim.sortentries();
                    auto stiffness = MatrixCSR( opl & ( vector_massmatrix & opr ) ); // MatrixCSR( stiffness_prelim );
                    
                    auto mass = MatrixCSR( incmatrix_t & ( scalar_massmatrix & incmatrix ) );
                    
                    auto stiffness_diagonal = SparseMatrix( DiagonalOperator( vector_massmatrix.diagonal() ) );
                    assert( stiffness_diagonal.issquare() );
                    assert( stiffness_diagonal.getdimin() == opr.getdimout() );
                    auto simplified_stiffness = MatrixCSR( opl & ( stiffness_diagonal & opr ) );
                    
                    auto idea_prelim = opl & opr;
                    idea_prelim.sortentries();
                    auto idea = MatrixCSR( idea_prelim );
                    
                    const auto& mat = stiffness;
//                     const auto& mat = simplified_stiffness;
                    
                    if(false)
                    {
                        cout << "average on the diagonal at level " << l << " is " << mat.diagonal().average() << nl;
                        contable << mat.diagonal().average() << nl;
                        contable.print( std::cout );
                    }
                    
                    assert( mat.diagonal().isfinite() );
                    
                    std::cout << mat << endl;
                    
                    {

                        FloatVector sol_original( opr.getdimin(), 0. );
                        FloatVector rhs_original( opr.getdimin(), 0. );
                        
                        sol_original.random(); sol_original.normalize( mass );
                        rhs_original.random(); rhs_original.normalize( mass );
                        
                        {
                            cout << "Filter out from x (CGM)" << endl;
                        
                            FloatVector sol( sol_original );
                            FloatVector rhs( rhs_original.getdimension(), 0. );
                            FloatVector residual( rhs );
                            
                            assert( sol.isfinite() );
                            
                            
                            ConjugateGradientSolverCSR( 
                                sol.getdimension(), 
                                sol.raw(), 
                                rhs.raw(), 
                                mat.getA(), mat.getC(), mat.getV(),
                                residual.raw(),
                                machine_epsilon,
                                1
                            );
                            sol.normalize( mass );
                            
                            assert( sol.isfinite() );
                            
                            std::cout << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm( mass ) << std::endl;
                            std::cout << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm( mass ) << std::endl;
                            
                            std::cout << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << std::endl;
                            std::cout << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << std::endl;
                            
                            contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
                            
                            
                            FloatVector sol2( sol_original );
                            sol2.random();
                            FloatVector rhs2( rhs_original.getdimension(), 0. );
                            FloatVector residual2( rhs );
                            
                            ConjugateGradientSolverCSR( 
                                sol2.getdimension(), 
                                sol2.raw(), 
                                rhs2.raw(), 
                                mat.getA(), mat.getC(), mat.getV(),
                                residual2.raw(),
                                machine_epsilon,
                                1
                            );
                            sol2.normalize( mass );
                            
                            
                            
                            
                            if( r == 1 ) {
                        
                                fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                    
                                VTKWriter vtk( M, fs, getbasename(__FILE__) );
                                vtk.writeCoordinateBlock();
                                vtk.writeTopDimensionalCells();
                                
                                vtk.writeVertexScalarData( sol,  "data1" , 1.0 );
                                vtk.writeVertexScalarData( sol2, "data2" , 1.0 );
                                // vtk.writeCellVectorData( interpol_grad, "gradient_interpolation" , 0.1 );
                                
                                fs.close();
                        
                            }

                            
                            
                        }

                        if(false)
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
                            sol.normalize( mass );
                            
                            std::cout << "\t\t\t x_0:       " << sol_original.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ax_0:      " << ( mat * sol_original ).norm( mass ) << std::endl;
                            std::cout << "\t\t\t b - Ax_0:  " << ( mat * sol_original - rhs ).norm( mass ) << std::endl;
                            
                            std::cout << "\t\t\t x:         " << sol.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ax:        " << ( mat * sol ).norm( mass ) << std::endl;
                            std::cout << "\t\t\t b - Ax:    " << ( mat * sol - rhs ).norm( mass ) << std::endl;
                            
                            contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
                        }

                        if(false)
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
                            residual.normalize( mass );
                            
                            std::cout << "\t\t\t b:       " << rhs_original.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ab:      " << ( mat * rhs ).norm( mass ) << std::endl;
                            
                            std::cout << "\t\t\t r:       " << residual.norm( mass ) << std::endl;
                            std::cout << "\t\t\t Ar:      " << ( mat * residual ).norm( mass ) << std::endl;
                            
                            std::cout << "\t\t\t Ar:      " << ( mat * ( rhs - mat * sol ) ).norm( mass ) << std::endl;
                            
                            contable << sol.norm( mass ) << ( mat * sol ).norm( mass );
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
