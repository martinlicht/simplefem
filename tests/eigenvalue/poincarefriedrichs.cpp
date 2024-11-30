

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../mesh/spanning.hpp"
#include "../../mesh/shelling.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: 3D curl estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = CrossedBricks3D();
    M.getcoordinates().scale( Constants::pi );
    
    M.check();
    
    // M.automatic_dirichlet_flags();
    // M.check_dirichlet_flags();

    
    // LOG << "Prepare scalar fields for testing..." << nl;
    

    if(false)
    {
        const auto shellings_found = generateShellings( M );

        LOG << shellings_found.size() << nl;

        for( int t = 0; t < shellings_found.size(); t++ )
        {
            const auto& shelling = shellings_found[t];
            LOG << t << "\t:\t";
            for( const int s : shelling ) LOG << s << space;
            LOG << nl;
        }

        return 0;
    }
            
    
    {
        LOG << "Generate spanning trees..." << nl;
        const std::pair< std::vector<int>, std::vector<std::vector<int>> > index2face_and_trees = list_face_spanning_trees( M );
        const auto& index2face = index2face_and_trees.first;
        const auto& trees      = index2face_and_trees.second;

        LOG << "Print spanning trees..." << nl;
        for( int t = 0; t < trees.size(); t++ )
        {
            const auto& tree = trees[t];
            LOG << t << "\t:\t";
            for( const int e : tree ) LOG << index2face[e] << space;
            LOG << nl;
        }

        return 0;

        
        LOG << "Generate shellings based on those trees..." << nl;
        
        for( const auto tree : trees )
        {
            assert( M.count_simplices(3) == tree.size()+1 );

            std::vector<std::vector<int>> volume_acceptable_face_list( M.count_simplices(3) );

            for( int index : tree ) 
            {
                int face = index2face[ index ];

                assert( 0 <= face and face < M.count_simplices(2) );

                const auto& parents = M.getsupersimplices(3,2,face);

                assert( parents.size() == 2 );

                for( int p : parents ) {
                    assert( 0 <= p and p < M.count_simplices(3) );
                    volume_acceptable_face_list[p].push_back(face);
                }
            
            }

            for( auto fl : volume_acceptable_face_list )
            {
                for( auto f : fl ) LOG << f << space;
                LOG << ", ";
            }
            LOG << nl;

            

            
            const auto shellings_found = generateShellings( M, volume_acceptable_face_list );

            LOG << shellings_found.size() << nl;

            for( int t = 0; t < shellings_found.size(); t++ )
            {
                const auto& shelling = shellings_found[t];
                LOG << t << "\t:\t";
                for( const int s : shelling ) LOG << s << space;
                LOG << nl;
            }

        }
        




        // LOG << "Print ordered spanning trees..." << nl;
        // for( int t = 0; t < trees.size(); t++ )
        // {
        //     LOG << t << "\t";
        //     const auto& tree = trees[t];
        //     const auto& results = list_ordered_face_spanning_trees( M, index2face, trees[t] );
        //     LOG << results.size() << nl;
        // }
        
    }


    return 0;

    LOG << M.text() << nl;

    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of curl operator (Sullivan)" << nl;

    const int min_l = 0; 
    const int max_l = 2;
    
    const int min_r = 1;
    const int max_r = 3;
    
    
    std::vector<ConvergenceTable> contables(max_r-min_r+1); //();
    for( int r = min_r; r <= max_r; r++ ){
        contables[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables[r-min_r] << "eigenvalue" << "ratio" << "log_2(ratio)" << "diff" << "log_2(diff)" << "u_mass" << "du_mass" << "time" << nl;
    } 

    

    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
      
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        if( l != 0 )
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "... assemble mass matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix;
            mat_A.sortandcompressentries();
                
            auto mat_Bt = vector_incmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
            mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & pseudo_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
            mat_B.sortandcompressentries();
            
            LOG << "... compose CSR system matrices" << nl;
    
            auto A  = MatrixCSR( mat_A  );
            auto Bt = MatrixCSR( mat_Bt );
            auto B  = MatrixCSR( mat_B  );
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // TODO: develop preconditioners 
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );

            auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                      + 
                      MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & pseudo_massmatrix & vector_diffmatrix & vector_incmatrix );
            auto PC = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )
                      + 
                      MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix );
            // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        

            LOG << "...begin inverse iteration" << nl;
            
            const int max_attempts = 1;

            for( int s = 0; s < max_attempts; s++ )
            {

                FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                candidate.random(); 
                candidate = A * candidate;
                candidate.normalize( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                
                const int max_inverseiterations = 10;

                Float newratio = -1;
                
                timestamp start = timestampnow();

                for( int t = 0; t < max_inverseiterations; t++ )
                {
                    
                    
                    // find the next candidate

                    FloatVector sol( A.getdimout(), 0. );
                    FloatVector aux( B.getdimout(), 0. );

                    const FloatVector rhs_sol = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate;
                    const FloatVector rhs_aux = FloatVector( B.getdimout(), 0. );

                    // auto residual = sol;
                    // ConjugateResidualSolverCSR_textbook( 
                    //     sol.getdimension(), 
                    //     sol.raw(), 
                    //     rhs_sol.raw(), 
                    //     A.getA(), A.getC(), A.getV(),
                    //     residual.raw(),
                    //     desired_precision,
                    //     -1
                    // );
                    
                    // const auto PAinv = pinv(PA,desired_precision,-1);
                    // const auto PCinv = pinv(PC,desired_precision,-1);
                    // sol = PAinv * rhs_sol;
                    
                    const auto PAinv = inv(PA,desired_precision,-1);
                    const auto PCinv = inv(PC,desired_precision,-1);
                    BlockHerzogSoodhalterMethod( 
                        sol, 
                        aux, 
                        rhs_sol, 
                        rhs_aux, 
                        A, Bt, B, C, 
                        desired_precision * sqrt(desired_precision),
                        -1,
                        // IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                        PAinv, PCinv
                    );

                    candidate = sol;
                    
                    
                    // assess the current candidate 

                    const auto A_candidate = A * candidate;
                    const auto M_candidate = ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ) * candidate; 
                    
                    const auto candidate_A_product = candidate * A_candidate; 
                    const auto candidate_M_product = candidate * M_candidate; 

                    newratio = candidate_A_product / candidate_M_product;

                    candidate /= sqrt(candidate_M_product); // Optional step

                    LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                    Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    
                    LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                    
                }

                timestamp end = timestampnow();

                // ... computed the solution

                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                
                LOG << "...compute error and residual" << nl;

                auto sol = candidate; 

                Float u_massnorm     = sol * ( vector_incmatrix_t * ( vector_massmatrix * vector_incmatrix * sol  ) );
                Float ucurl_massnorm = sol * ( mat_A * sol );
                Float curratio       = ucurl_massnorm / u_massnorm;
                Float u_defectmass   = ( B * sol ).norm_sq( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                
                LOG << "ratio:           " << newratio << nl;
                LOG << "ratio:           " << curratio << nl;
                LOG << "u mass:          " << u_massnorm << nl;
                LOG << "u curl mass      " << ucurl_massnorm << nl;
                LOG << "u defect mass:   " << u_defectmass << nl;
                
                const Float true_eigenvalue = 2.; // 3.0 is the true value 

                contables[r-min_r] << newratio;
                contables[r-min_r] << newratio / true_eigenvalue - 1.;
                contables[r-min_r] << - std::log2( newratio / true_eigenvalue - 1. ); 
                contables[r-min_r] << newratio - true_eigenvalue;
                contables[r-min_r] << std::log2( newratio - true_eigenvalue );
                contables[r-min_r] << u_massnorm;
                contables[r-min_r] << ucurl_massnorm;
                contables[r-min_r] << Float( end - start );
                contables[r-min_r] << nl;

                contables[r-min_r].lg();
            
            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        

    } 

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
