

/**/

#include <cmath>

#include <algorithm>
#include <fstream>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"



Float alpha(const Float x) { 
    Float maxnorm = absolute(x);
    return std::exp( 1 - 1. / maxnorm );
};

Float alpha_dev(const Float x) { 
    Float maxnorm = absolute(x);
    Float maxnorm_sq = square(maxnorm);            
    return std::exp( 1 - 1. / maxnorm ) / ( maxnorm_sq );
};

const Float B = 2.;

FloatVector trafo(const FloatVector& vec) { 
    assert( vec.getdimension() == 3 );
    
    return vec / B + alpha( vec.maxnorm() ) * ( 1./vec.l2norm() - 1./vec.maxnorm()/B ) * vec;
};

        
DenseMatrix jacobian(const FloatVector& vec) { 
    assert( vec.getdimension() == 3 );
    
    // return DenseMatrix( 3, 3, kronecker<int> );

    Float x = vec[0];
    Float y = vec[1];
    Float z = vec[2];

    Float ax = absolute(x);
    Float ay = absolute(y);
    Float az = absolute(z);

    Float sx = sign(x);
    Float sy = sign(y);
    Float sz = sign(z);

    Float dx = ( ax >= ay and ax >= az ) ? sx : 0.;
    Float dy = ( ay >= ax and ay >= az ) ? sy : 0.;
    Float dz = ( az >= ax and az >= ay ) ? sz : 0.;
    
    Float li  = maximum(ax,maximum(ay,az));
    Float lis = li*li;
    Float l2s = x*x + y*y + z*z;
    Float l2  = sqrt( l2s );
    Float l2c = l2*l2s;

    Float a     = std::exp( 1. - 1. / li );
    Float a_dev = a / lis;
    
    return DenseMatrix( 3, 3, {
        1./B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sx * dx * ( 1/l2 - 1/li/B ) + a * ( -x/l2c - sx*dx/lis/B ) ) * x 
        ,
                                       (                                       a * ( -y/l2c               ) ) * x 
        ,
                                       (                                       a * ( -z/l2c               ) ) * x 
        ,
        //
                                       (                                       a * ( -x/l2c               ) ) * y 
        ,
        1./B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sy * dy * ( 1/l2 - 1/li/B ) + a * ( -y/l2c - sy*dy/lis/B ) ) * y 
        ,
                                       (                                       a * ( -z/l2c               ) ) * y 
        ,
        //
                                       (                                       a * ( -x/l2c               ) ) * z 
        ,
                                       (                                       a * ( -y/l2c               ) ) * z 
        ,
        1./B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sz * dz * ( 1/l2 - 1/li/B ) + a * ( -z/l2c - sz*dz/lis/B ) ) * z 
        ,
    });
};

        
FloatVector physical_f(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    return FloatVector({ 
        1.
        });
};

FloatVector parametric_f(const FloatVector& vec) 
{
    return physical_f( trafo(vec) );
};


DenseMatrix physical_A(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    Float scalar = vec.l2norm() > 0.5 ? 1. : 2. ;
    return scalar * DenseMatrix(3,3, kronecker<int> );
};

DenseMatrix parametric_A(const FloatVector& vec) 
{
    return physical_A( trafo(vec) );
};





DenseMatrix weight_scalar(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    // return DenseMatrix(1,1, kronecker<int> );
    return DenseMatrix(1,1,absolute(Determinant(jacobian(vec))));
};
        
DenseMatrix weight_vector(const FloatVector& vec)
{
    assert( vec.getdimension() == 3 );
    // return DenseMatrix(3,3, kronecker<int> );
    auto jac = jacobian(vec);

    auto det = Determinant(jac);

    auto invjac = Inverse( jac );

    return absolute(det) * invjac * parametric_A(vec) * Transpose(invjac);
};
        
        



FloatVector IncreaseResolution( const MeshSimplicial3D& mesh, const FloatVector& lores )
{
    int counter_vertices = mesh.count_vertices();
    int counter_edges    = mesh.count_edges();

    assert( lores.getdimension() == counter_vertices );

    FloatVector hires( counter_vertices + counter_edges );

    for( int v = 0; v < counter_vertices; v++ ) hires[v] = lores[v];

    for( int e = 0; e < counter_edges; e++ ) {
        int v0 = mesh.get_edge_vertex(e,0);
        int v1 = mesh.get_edge_vertex(e,1);
        Float value = ( lores[v0] + lores[v1] ) / 2.;
        hires[counter_vertices + e] = value;
    }

    return hires;
}




using namespace std;

int main()
{
        
    LOG << "Unit Test for transformed 3D Poisson Problem" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = FicheraCorner3D();
    
    M.check();

    
    for( int e = 0; e < M.count_edges(); e++ ) {
        Float D = M.getDiameter(1,e);
        if( 1.9 < D and D < 2.1 )
            M.bisect_edge(e); 
    }

    
    for( int s = 0; s < M.count_simplices(2); s++ ) 
    {
        if( M.getsupersimplices(3,2,s).size() > 1 ) continue;
        
        auto midpoint = M.get_midpoint(2,s);

        bool eligible = ( midpoint[0] < 0. );

        if( not eligible ) continue;

        M.set_flag( 2, s, SimplexFlagDirichlet );
        
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 0 ), SimplexFlagDirichlet );
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 1 ), SimplexFlagDirichlet );
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 2 ), SimplexFlagDirichlet );
        
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 0 ), SimplexFlagDirichlet );
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 1 ), SimplexFlagDirichlet );
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 2 ), SimplexFlagDirichlet );
        
    }
    
    M.check_dirichlet_flags(false);

    // for( int s = 0; s < M.count_simplices(2); s++ ) 
    // {
    //     auto mid = M.get_midpoint(2,s);
    //     bool b = ( mid[0] < 0. );
    //     if( M.get_flag(d,s) == SimplexFlagDirichlet and b )
    //         M.set_flag( d, s, SimplexFlagNull );
    // }
        
    M.check_dirichlet_flags(false);

    
    LOG << "Prepare scalar fields for testing..." << nl;


    

    
    
    
    

    

    LOG << "Solving Poisson Problem with Neumann boundary conditions" << nl;

    const int min_l = 1; 
    const int max_l = 3;

    const int min_r = 1;
    const int max_r = 1;

    const int r_plus = 1;
    const int w_plus = 1;
    
    ConvergenceTable contable("Mass error");
    
    contable << "u_error" << "du_error" << "residual" << "time" << nl;
    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
    
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        // if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            int w = r;

            LOGPRINTF("Polynomial degrees: r=%d w=%d r+=%d w+=%d \n", r, w, r_plus, w_plus );
            
            LOG << "...assemble scalar mass matrices" << nl;
    
            auto     scalar_massmatrix = MatrixCSR( FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 0, r         , w         , weight_scalar ) );
            auto aug_scalar_massmatrix = MatrixCSR( FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 0, r + r_plus, w + w_plus, weight_scalar ) );
            // auto     scalar_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r          ) );
            // auto aug_scalar_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus ) );

            LOG << "...assemble vector mass matrix" << nl;
    
            auto     vector_massmatrix = MatrixCSR( FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 1, r-1         , w         , weight_vector ) );
            auto aug_vector_massmatrix = MatrixCSR( FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 1, r-1 + r_plus, w + w_plus, weight_vector ) );
            // auto     vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1          ) );
            // auto aug_vector_massmatrix = MatrixCSR( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 + r_plus ) );
            
            LOG << "...assemble differential matrix and transpose" << nl;

            auto     diffmatrix = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r          ) );
            auto aug_diffmatrix = MatrixCSR( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r + r_plus ) );

            auto     diffmatrix_t =     diffmatrix.getTranspose();
            auto aug_diffmatrix_t = aug_diffmatrix.getTranspose();

            LOG << "...assemble inclusion matrix and transpose" << nl;
    
            auto     incmatrix = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r          ) );
            auto aug_incmatrix = MatrixCSR( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r + r_plus ) );
            
            auto     incmatrix_t =     incmatrix.getTranspose();
            auto aug_incmatrix_t = aug_incmatrix.getTranspose();

            LOG << "...assemble elevation matrix" << nl;
    
            auto elevation_matrix = MatrixCSR( FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus ) );

            display_mallinfo();

            {

                const auto& function_rhs  = parametric_f;
                
                LOG << "...interpolate explicit solution and rhs" << nl;
    
                FloatVector     interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,        function_rhs  );
                FloatVector aug_interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r+r_plus, function_rhs  );

                FloatVector     rhs =     incmatrix_t * (     scalar_massmatrix *     interpol_rhs );
                FloatVector aug_rhs = aug_incmatrix_t * ( aug_scalar_massmatrix * aug_interpol_rhs );

                FloatVector     sol(     incmatrix.getdimin(), 0. );
                FloatVector aug_sol( aug_incmatrix.getdimin(), 0. );
                
                FloatVector     residual(     incmatrix.getdimin(), 0. );
                FloatVector aug_residual( aug_incmatrix.getdimin(), 0. );
                
                timestamp solver_time;

                LOG << "...iterative solver" << nl;
                
                // if(false)
                {
                    LOG << "[0]" << nl;
                    auto opr  = (diffmatrix) & (incmatrix);
                    LOG << "[1]" << nl;
                    // auto opl  = MatrixCSR(incmatrix_t) & MatrixCSR(diffmatrix_t);
                    auto opl  = opr.getTranspose();
                    LOG << "[2]" << nl;
                    auto stiffness = (opl) & (vector_massmatrix) & (opr);
                    LOG << "[3]" << nl;
                    auto& stiffness_csr = stiffness; //MatrixCSR( stiffness );
                    LOG << "[4]" << nl;
                    
                    sol.zero();

                    auto diagonal = stiffness.diagonal();
                    // FloatVector residual( rhs );
                
                    timestamp start = gettimestamp();
                    ConjugateGradientSolverCSR_SSOR( 
                        sol.getdimension(), 
                        sol.raw(), 
                        rhs.raw(), 
                        stiffness_csr.getA(), stiffness_csr.getC(), stiffness_csr.getV(),
                        residual.raw(),
                        desired_precision,
                        0,
                        diagonal.raw(),
                        1.0
                    );
                    timestamp end = gettimestamp();
                    solver_time = end - start;
                    LOG << "\t\t\t Time: " << timestamp2measurement( solver_time ) << nl;

                    residual = rhs - stiffness * sol;

                    display_mallinfo();

                }

                {
                    LOG << "[0]" << nl;
                    auto aug_opr  = (aug_diffmatrix) & (aug_incmatrix);
                    LOG << "[1]" << nl;
                    // auto aug_opl  = MatrixCSR(aug_incmatrix_t) & MatrixCSR(aug_diffmatrix_t);
                    auto aug_opl  = aug_opr.getTranspose();
                    LOG << "[2]" << nl;
                    auto aug_stiffness = (aug_opl) & (aug_vector_massmatrix) & (aug_opr);
                    LOG << "[3]" << nl;
                    auto& aug_stiffness_csr = aug_stiffness; //MatrixCSR( aug_stiffness );
                    LOG << "[4]" << nl;

                    aug_sol.zero();

                    auto aug_diagonal = aug_stiffness.diagonal();
                    // FloatVector aug_residual( aug_rhs );
                    
                    timestamp start = gettimestamp();
                    ConjugateGradientSolverCSR_SSOR( 
                        aug_sol.getdimension(), 
                        aug_sol.raw(), 
                        aug_rhs.raw(), 
                        aug_stiffness_csr.getA(), aug_stiffness_csr.getC(), aug_stiffness_csr.getV(),
                        aug_residual.raw(),
                        desired_precision,
                        0,
                        aug_diagonal.raw(),
                        1.0
                    );
                    timestamp end = gettimestamp();
                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                    aug_residual = aug_rhs - aug_stiffness * aug_sol;

                    display_mallinfo();

                }

            
                
                
                LOG << "...compute error and residual:" << nl;

                FloatVector error     = aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol;
                FloatVector graderror = aug_diffmatrix * ( aug_incmatrix * aug_sol - elevation_matrix * incmatrix * sol );
                Float errornorm       = std::sqrt( error * ( aug_scalar_massmatrix * error ) );
                Float graderrornorm   = std::sqrt( graderror * ( aug_vector_massmatrix * graderror ) );
                Float residualnorm    = residual.norm();

                LOG << "error:     " << errornorm    << nl;
                LOG << "graderror: " << graderrornorm << nl;
                LOG << "residual:  " << residualnorm << nl;
                LOG << "time:      " << Float( solver_time ) << nl;
                        
                contable << errornorm;
                contable << graderrornorm;
                contable << residualnorm;
                contable << Float( solver_time );
                contable << nl;
                    
                contable.lg();


                if( r == 1 ){
                    
                    auto computed_grad = diffmatrix * incmatrix * sol;
                    
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    vtk.writeCoordinateBlock();
                    vtk.writeTopDimensionalCells();
                    vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    vtk.writeCellVectorData( computed_grad, "gradient_interpolation" , 0.1 );
                    fs.close();
                }


            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        
        

    } 
    
    
    LOG << "Finished Unit Test" << nl;
    
    return 0;
}
