

/**/

#include <cmath>

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../mesh/neumannestimate.hpp"
#include "../../mesh/spanning.hpp"
#include "../../mesh/shelling.hpp"
#include "../../mesh/shelling2.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.avgsullivan.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/utilities.hpp"


// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D curl PF estimates using shellings" << nl;
    
    LOG << "Initial mesh..." << nl;

    const bool do_gradient   = true;
    const bool do_curl       = true;
    const bool do_divergence = true;
    
    MeshSimplicial3D M = UnitCube3D();
    // MeshSimplicial3D M = StandardCubeFive3D();
    // MeshSimplicial3D M = CrossedBricks_Five3D();
    // MeshSimplicial3D M = CrossedBricks3D();
    // MeshSimplicial3D M = FicheraCorner3D();
    // MeshSimplicial3D M = RandomPolyhedralSphere(0);
    M.check();
    M.getCoordinates().shake_random();

    // return 0;

    LOG << "Number of tetrahedra: " << M.count_tetrahedra() << nl;
    
    Float PF_estimate_via_shellings[3] = { notanumber, notanumber, notanumber };

    Float grad_estimate = NeumannEstimate( M );

    LOG << grad_estimate << nl;

    
    
    // if(false)
    for( int k = 0; k < 3; k++ )
    {
        
        PF_estimate_via_shellings[k] = std::numeric_limits<Float>::infinity();

        auto shellings_found = generate_shellings2( M, k );

        LOG << shellings_found.size() << nl;

        // typedef decltype(shellings_found[0]) shelling;

        std::sort( shellings_found.begin(), shellings_found.end(), 
                    []( const shelling& s1, const shelling& s2 ) -> bool{ return s1.weight_reflection < s2.weight_reflection; } 
                    );

        for( int t = 0; t < shellings_found.size(); t++ )
        {
            const auto& shelling = shellings_found[t];
            LOG << t << "\t:\t";
            for( const auto& s : shelling ) { LOG << s << space; }
            LOG << "\t" << shelling.weight_reflection << nl;
        }

        PF_estimate_via_shellings[k] = shellings_found.front().weight_reflection;
            
    }

    return 0;
    
            
    /*
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
            
    
    if(false)
    {
        const auto shellings_found = generate_ranked_shelling( M );

        LOG << shellings_found.size() << nl;

        for( int t = 0; t < shellings_found.size(); t++ )
        {
            Float weight = 0.;
            const auto& shelling = shellings_found[t];
            LOG << t << "\t:\t";
            for( const auto& s : shelling ) { LOG << s.first << " (" << s.second << ")" << space; weight += s.second; }
            LOG << space << weight;
            LOG << nl;

            std::vector<int> indices( shelling.size() );
            for( int i = 0; i < shelling.size(); i++ ) indices[i] = shelling[i].first;
            Float estimate = estimate_shelling_quality( M, indices, 1 );

            PF_estimate_via_shellings[0] = minimum( PF_estimate_via_shellings[0], estimate );
        }

        LOG << "PF estimate via shellings: " << PF_estimate_via_shellings << nl;

        return 0;
    }
    */
            
    
    if(false)
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
        
        for( const auto& tree : trees )
        {
            assert( M.count_simplices(3) == tree.size()+1 );

            std::vector<bool> acceptable_face_list( M.count_simplices(2) );

            for( int index : tree ) 
            {
                int face = index2face[ index ];

                assert( 0 <= face and face < M.count_simplices(2) );

                acceptable_face_list[face] = true;
            
            }

            for( auto fl : acceptable_face_list ) LOG << (bool)fl << space;
            LOG << nl;

            

            
            const auto shellings_found = generateShellings( M, acceptable_face_list );

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


    // LOG << M.text() << nl;
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of curl operator (Whitney)" << nl;

    const int min_l = 3; 
    const int max_l = 3;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
    std::vector<ConvergenceTable> contables_scalar(max_r-min_r+1);
    std::vector<ConvergenceTable> contables_vector(max_r-min_r+1);
    std::vector<ConvergenceTable> contables_pseudo(max_r-min_r+1);
    
    for( int r = min_r; r <= max_r; r++ )
    {
        contables_scalar[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables_scalar[r-min_r] << "eigenvalue (iterated)" << "eigenvalue (post)" << "PF (post)" << "PF ratio" << "u_mass" << "du_mass" << "time" << nl;

        contables_vector[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables_vector[r-min_r] << "eigenvalue (iterated)" << "eigenvalue (post)" << "PF (post)" << "PF ratio" << "u_mass" << "du_mass" << "time" << nl;

        contables_pseudo[r-min_r].table_name = "Mass error and numerical residuals r=" + std::to_string(r);
        contables_pseudo[r-min_r] << "eigenvalue (iterated)" << "eigenvalue (post)" << "PF (post)" << "PF ratio" << "u_mass" << "du_mass" << "time" << nl;
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
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix pseudo_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r   );
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 );

            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r   );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix pseudo_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r   );
            SparseMatrix pseudo_incmatrix_t = pseudo_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            SparseMatrix pseudo_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r );
            SparseMatrix pseudo_diffmatrix_t = pseudo_diffmatrix.getTranspose();


            SparseMatrix scalar_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r-1, 1 );
            SparseMatrix scalar_elevationmatrix_t = scalar_elevationmatrix.getTranspose();

            SparseMatrix vector_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1 );
            SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

            SparseMatrix pseudo_elevationmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r-1, 1 );
            SparseMatrix pseudo_elevationmatrix_t = pseudo_elevationmatrix.getTranspose();


            LOG << "averaging matrix for scalar fields" << nl;

            SparseMatrix scalar_averaging = FEECSullivanAveragingMatrix( M, M.getinnerdimension(), 0, r, FEECAveragingMode::weighted_uniformly );

            




            if( do_gradient ) {
                
                LOG << "... compose system matrices (SCALAR)" << nl;

                // PING;
                // scalar_diffmatrix & scalar_incmatrix;
                // PING;
                // vector_elevationmatrix & scalar_diffmatrix;
                // PING;
                // vector_massmatrix & scalar_elevationmatrix;
                // PING;
                auto temp_scalar = vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix;
                auto temp_scalar_t = temp_scalar.getTranspose();
                auto mat_A  = temp_scalar_t & vector_massmatrix & temp_scalar;
                mat_A.sortandcompressentries();
                    
                LOG << "... compose CSR system matrices (SCALAR)" << nl;
        
                auto A  = MatrixCSR( mat_A  );
                
                            

                LOG << "... begin inverse iteration" << nl;
                
                const int max_attempts = 1;

                for( int s = 0; s < max_attempts; s++ )
                {

                    FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                    candidate.random(); 
                    candidate = A * candidate;
                    candidate.normalize( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    
                    const int max_inverseiterations = 10;

                    Float newratio = -1;
                    
                    FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, 
                    [](const FloatVector& vec) -> FloatVector { return FloatVector({ 1. }); }
                    );

                    FloatVector conforming_one = scalar_averaging * interpol_one;

                    Float mass_of_conforming_one  = ( scalar_massmatrix * interpol_one ) * interpol_one;

                    timestamp start = timestampnow();
                    
                    for( int t = 0; t <= max_inverseiterations; t++ )
                    {
                        
                        
                        // find the next candidate

                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * candidate ) * interpol_one;
                            candidate = candidate - ( average / mass_of_conforming_one ) * conforming_one;
                        }

                        FloatVector sol( A.getdimout(), 0. );
                        
                        FloatVector rhs_sol = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate;
                        
                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * rhs_sol ) * interpol_one;
                            rhs_sol = rhs_sol - ( average / mass_of_conforming_one ) * conforming_one;
                        }

                        auto residual = sol;
                        
                        ConjugateResidualSolverCSR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs_sol.raw(), 
                            A.getA(), A.getC(), A.getV(),
                            residual.raw(),
                            desired_precision,
                            -1
                        );

                        {
                            DenseMatrix Bt( A.getdimout(), 1, 1. );
                            DenseMatrix B = Transpose(Bt);
                            DenseMatrix C(1,1,0.);
                            
                            auto aux     = FloatVector(1,0.);
                            auto rhs_aux = FloatVector(1,0.);
                            
                            BlockHerzogSoodhalterMethod( 
                                sol, 
                                aux, 
                                rhs_sol, 
                                rhs_aux, 
                                A, Bt, B, C, 
                                desired_precision * std::sqrt(desired_precision),
                                -1,
                                IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                            );

                        }

                        {
                            Float average = ( scalar_massmatrix * scalar_incmatrix * sol ) * interpol_one;
                            sol = sol - ( average / mass_of_conforming_one ) * conforming_one;
                        }
                        
                        candidate = sol;
                        
                        
                        // assess this new candidate 

                        const auto A_candidate = A * candidate;
                        const auto M_candidate = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate; 
                        
                        const auto candidate_A_product = candidate * A_candidate; 
                        const auto candidate_M_product = candidate * M_candidate; 

                        newratio = candidate_A_product / candidate_M_product;

                        candidate /= std::sqrt(candidate_M_product); // Optional step

                        LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;
                        
                        Float u_residualmass_sq   = ( A * sol - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                        
                        LOG << "current residual: " << u_residualmass_sq << nl;

                        
                    }

                    timestamp end = timestampnow();

                    // ... computed the solution

                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    LOG << "... compute error and residual" << nl;

                    auto sol = candidate; 

                    Float u_massnorm     = sol * ( scalar_incmatrix_t * ( scalar_massmatrix * scalar_incmatrix * sol  ) );
                    Float ugrad_massnorm = sol * ( mat_A * sol );
                    Float curratio       = ugrad_massnorm / u_massnorm;
                    
                    LOG << "ratio:       " << newratio << nl;
                    LOG << "ratio:       " << curratio << nl;
                    LOG << "u mass:      " << u_massnorm << nl;
                    LOG << "u grad mass: " << ugrad_massnorm << nl;

                    LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                    
                    const Float true_eigenvalue = 1.;
                    // 1.0 is the true value 
                    // 3.0 is the true value 

                    contables_scalar[r-min_r] << newratio;
                    contables_scalar[r-min_r] << curratio;
                    contables_scalar[r-min_r] << std::sqrt(curratio);
                    contables_scalar[r-min_r] << PF_estimate_via_shellings[1] * std::sqrt(curratio);
                    contables_scalar[r-min_r] << u_massnorm;
                    contables_scalar[r-min_r] << ugrad_massnorm;
                    contables_scalar[r-min_r] << Float( end - start );
                    contables_scalar[r-min_r] << nl;

                    contables_scalar[r-min_r].lg();
                
                }
                
            }







            if( do_curl ) {
                
                LOG << "... compose system matrices (VECTOR)" << nl;
    
                auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & pseudo_elevationmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix;
                mat_A.sortandcompressentries();
                    
                auto mat_Bt = vector_incmatrix_t & vector_massmatrix & vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); // lower left
                mat_B.sortandcompressentries();
                
                LOG << "... compose CSR system matrices (VECTOR)" << nl;
        
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                // auto PA = IdentityMatrix( A.getdimin() );
                // auto PC = IdentityMatrix( C.getdimin() );

                auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                          + 
                          MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & pseudo_elevationmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix );
                
                auto PC = MatrixCSR( scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix )
                          + 
                          MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_elevationmatrix & scalar_diffmatrix & scalar_incmatrix );
                // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                            

                LOG << "... begin inverse iteration (VECTOR)" << nl;
            
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
                            desired_precision * std::sqrt(desired_precision),
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

                        candidate /= std::sqrt(candidate_M_product); // Optional step

                        LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                        Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                        Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                        
                        LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                        
                    }

                    timestamp end = timestampnow();

                    // ... computed the solution

                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    LOG << "... compute error and residual (VECTOR)" << nl;

                    auto sol = candidate; 

                    Float u_massnorm     = sol * ( vector_incmatrix_t * ( vector_massmatrix * vector_incmatrix * sol  ) );
                    Float ucurl_massnorm = sol * ( mat_A * sol );
                    Float curratio       = ucurl_massnorm / u_massnorm;
                    Float u_defectmass   = ( B * sol ).norm_sq( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ); 
                    
                    LOG << "ratio (iteration):   " << newratio << nl;
                    LOG << "ratio (postprocess): " << curratio << nl;
                    LOG << "u mass:              " << u_massnorm << nl;
                    LOG << "u curl mass          " << ucurl_massnorm << nl;
                    LOG << "u defect mass:       " << u_defectmass << nl;
                    
                    LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                        
                    contables_vector[r-min_r] << newratio;
                    contables_vector[r-min_r] << curratio;
                    contables_vector[r-min_r] << std::sqrt(curratio);
                    contables_vector[r-min_r] << PF_estimate_via_shellings[1] * std::sqrt(curratio);
                    contables_vector[r-min_r] << u_massnorm;
                    contables_vector[r-min_r] << ucurl_massnorm;
                    contables_vector[r-min_r] << Float( end - start );
                    contables_vector[r-min_r] << nl;

                    contables_vector[r-min_r].lg();
                
                }
            
            }







            if( do_divergence ) {
                
                LOG << "... compose system matrices (PSEUDO)" << nl;
    
                auto mat_A  = pseudo_incmatrix_t & pseudo_diffmatrix_t & volume_massmatrix & pseudo_diffmatrix & pseudo_incmatrix;
                mat_A.sortandcompressentries();
                    
                auto mat_Bt = pseudo_incmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix; // upper right
                mat_Bt.sortandcompressentries();
                
                auto mat_B = mat_Bt.getTranspose(); // lower left
                mat_B.sortandcompressentries();
                
                LOG << "... compose CSR system matrices (PSEUDO)" << nl;
        
                auto A  = MatrixCSR( mat_A  );
                auto Bt = MatrixCSR( mat_Bt );
                auto B  = MatrixCSR( mat_B  );
                
                auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                
                // auto PA = IdentityMatrix( A.getdimin() );
                // auto PC = IdentityMatrix( C.getdimin() );

                auto PA = MatrixCSR( pseudo_incmatrix_t & pseudo_massmatrix & pseudo_incmatrix )
                          + 
                          MatrixCSR( pseudo_incmatrix_t & pseudo_diffmatrix_t & volume_massmatrix & pseudo_diffmatrix & pseudo_incmatrix );

                auto PC = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                          + 
                          MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & pseudo_elevationmatrix_t & pseudo_massmatrix & pseudo_elevationmatrix & vector_diffmatrix & vector_incmatrix );

                // LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
                // LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                            

                LOG << "... begin inverse iteration (PSEUDO)" << nl;
            
                const int max_attempts = 1;

                for( int s = 0; s < max_attempts; s++ )
                {

                    FloatVector candidate = FloatVector( A.getdimout(), 0. ); 
                    candidate.random(); 
                    candidate = A * candidate;
                    candidate.normalize( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ); 
                    
                    const int max_inverseiterations = 10;

                    Float newratio = -1;
                    
                    timestamp start = timestampnow();

                    for( int t = 0; t < max_inverseiterations; t++ )
                    {
                        
                        
                        // find the next candidate

                        FloatVector sol( A.getdimout(), 0. );
                        FloatVector aux( B.getdimout(), 0. );

                        const FloatVector rhs_sol = ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ) * candidate;
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
                            desired_precision * std::sqrt(desired_precision),
                            -1,
                            // IdentityMatrix( A.getdimin() ), IdentityMatrix( C.getdimin() ) 
                            PAinv, PCinv
                        );

                        candidate = sol;
                        
                        
                        // assess the current candidate 

                        const auto A_candidate = A * candidate;
                        const auto M_candidate = ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ) * candidate; 
                        
                        const auto candidate_A_product = candidate * A_candidate; 
                        const auto candidate_M_product = candidate * M_candidate; 

                        newratio = candidate_A_product / candidate_M_product;

                        candidate /= std::sqrt(candidate_M_product); // Optional step

                        LOG << "current ratio: " << newratio << " (" << t << "/" << max_inverseiterations << ")" << nl;

                        Float u_residualmass_sq   = ( A * sol + Bt * aux - rhs_sol ).norm_sq(); // ( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                        Float aux_residualmass_sq = ( B * sol            - rhs_aux ).norm_sq(); // ( pseudo_incmatrix_t * pseudo_massmatrix * pseudo_incmatrix ); 
                        
                        LOG << "current residuals: " << u_residualmass_sq << tab << aux_residualmass_sq << nl;

                        
                    }

                    timestamp end = timestampnow();

                    // ... computed the solution

                    LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                    
                    LOG << "... compute error and residual (PSEUDO)" << nl;

                    auto sol = candidate; 

                    Float u_massnorm     = sol * ( pseudo_incmatrix_t * ( pseudo_massmatrix * pseudo_incmatrix * sol  ) );
                    Float ucurl_massnorm = sol * ( mat_A * sol );
                    Float curratio       = ucurl_massnorm / u_massnorm;
                    Float u_defectmass   = ( B * sol ).norm_sq( vector_incmatrix_t * vector_massmatrix * vector_incmatrix ); 
                    
                    LOG << "ratio (iteration):   " << newratio << nl;
                    LOG << "ratio (postprocess): " << curratio << nl;
                    LOG << "u mass:              " << u_massnorm << nl;
                    LOG << "u curl mass          " << ucurl_massnorm << nl;
                    LOG << "u defect mass:       " << u_defectmass << nl;
                    
                    LOG << "PF constant estimates: " << 1./std::sqrt(curratio) << space  << 1./std::sqrt(newratio) << nl;
                        
                    contables_pseudo[r-min_r] << newratio;
                    contables_pseudo[r-min_r] << curratio;
                    contables_pseudo[r-min_r] << std::sqrt(curratio);
                    contables_pseudo[r-min_r] << PF_estimate_via_shellings[2] * std::sqrt(curratio);
                    contables_pseudo[r-min_r] << u_massnorm;
                    contables_pseudo[r-min_r] << ucurl_massnorm;
                    contables_pseudo[r-min_r] << Float( end - start );
                    contables_pseudo[r-min_r] << nl;

                    contables_pseudo[r-min_r].lg();
                
                }
            
            }


        }


        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        

    } 

    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}


