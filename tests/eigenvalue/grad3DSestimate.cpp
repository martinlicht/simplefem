

/**/

#include <cmath>
#include <cstdlib>

#include <limits>
#include <string>
#include <vector>

#include "../../base/include.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../fem/global.avgsullivan.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


#include <queue>
#include <utility> // for std::pair
#include <numeric> // for accumulate



// TODO(martinlicht): can this be a local lambda?
// Comparator for the priority queue (min-heap, based on weights)
struct Compare {
    bool operator()(const std::pair<int, Float>& p1, const std::pair<int, Float>& p2) {
        // Use > to prioritize the pair with the smaller weight
        return p1.second > p2.second;
    }
};










// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: 3D grad estimate" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitSimplex3D(); 
    M.getCoordinates().scale( 3. );
    M.getCoordinates().setdata(3, 2, 1./3. );
    
    // MeshSimplicial3D M = UnitCube3D(); 
    // MeshSimplicial3D M = FicheraCorner3D();
    // MeshSimplicial3D M = CrossedBricks3D();
    
    M.check();
    
    bool do_neumann = false;
    // if( not do_neumann ){
    //     M.automatic_dirichlet_flags();
    //     M.check_dirichlet_flags();
    // }

    {
        LOG << "Fine-tuned boundary conditions" << nl;

        int first_bc_face = 0;
        
        if( argc > 1 )
        {
            const char* end = nullptr;
            bool has_overflown;
            int value = string_to_integer( argv[1], &end, 10, has_overflown );

            // Check if the entire argument was parsed and within int range
            if( *end != '\0' ) {

                LOG << "Error: The provided argument is not a valid integer:" << argv[1] << "\n";

            } else if( value != static_cast<int>(value) ) {

                LOG << "Error: The provided argument is out of 'int' range.\n";

            } else {
                
                int number = static_cast<int>(value);
                LOG << "Command-line argument provided: " << number << nl;

                first_bc_face = number;

            }
        } else {
            LOG << "Dirichlet faces start, per default, at " << first_bc_face << nl;
        }

        for( int f = first_bc_face; f <= 3; f++ ) {
            M.set_flag( 2, f, SimplexFlag::SimplexFlagDirichlet );
        }
        M.complete_dirichlet_flags_from_facets();
        M.check_dirichlet_flags(false);

        if( first_bc_face > 3 ) do_neumann = true;
    
    }

    for( int f = 0; f < M.count_faces(); f++ )
    {
        LOG << f << ": ";
        for( int i = 0; i <= 2; i++ ) LOG << M.get_subsimplex( 2, 0, f, i ) << space;
        LOG << " > " << (int)M.get_flag( 2, f ) << nl;
    }
    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    
            
    
    
    
    
    

    

    LOG << "Estimating Poincare-Friedrichs constant of grad operator (Sullivan)" << nl;

    const int min_l = 3; 
    const int max_l = 3;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
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
            
            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            LOG << "averaging matrix for scalar fields" << nl;

            SparseMatrix scalar_averaging = FEECSullivanAveragingMatrix( M, M.getinnerdimension(), 0, r+1, FEECAveragingMode::weighted_uniformly );

            LOG << "... compose system matrices" << nl;

            auto mat_A  = scalar_incmatrix_t & scalar_diffmatrix_t & ( vector_massmatrix ) & scalar_diffmatrix & scalar_incmatrix;
            mat_A.sortandcompressentries();
                
            LOG << "... compose CSR system matrices" << nl;
    
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
                
                FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r+1, 
                    [](const FloatVector& vec) -> FloatVector { return FloatVector({ 1. }); }
                    );

                FloatVector conforming_one = scalar_averaging * interpol_one;

                Float mass_of_conforming_one  = ( scalar_massmatrix * interpol_one ) * interpol_one;

                timestamp start = timestampnow();
                
                for( int t = 0; t <= max_inverseiterations; t++ )
                {
                    
                    
                    // find the next candidate

                    if( do_neumann )
                    {
                        Float average = ( scalar_massmatrix * scalar_incmatrix * candidate ) * interpol_one;
                        candidate = candidate - ( average / mass_of_conforming_one ) * conforming_one;
                    }

                    FloatVector sol( A.getdimout(), 0. );
                    
                    FloatVector rhs_sol = ( scalar_incmatrix_t * scalar_massmatrix * scalar_incmatrix ) * candidate;
                    
                    if( do_neumann )
                    {
                        Float average = ( scalar_massmatrix * scalar_incmatrix * rhs_sol ) * interpol_one;
                        rhs_sol = rhs_sol - ( average / mass_of_conforming_one ) * conforming_one;
                    }

                    auto residual = sol;
                    
                    if( not do_neumann ) {

                        ConjugateResidualSolverCSR( 
                            sol.getdimension(), 
                            sol.raw(), 
                            rhs_sol.raw(), 
                            A.getA(), A.getC(), A.getV(),
                            residual.raw(),
                            desired_precision,
                            -1
                        );

                    } else {
                        
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

                    if( do_neumann )
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
                
                const Float true_eigenvalue = ( do_neumann ? 1.0 : 3.0 );
                // 1.0 is the true value 
                // 3.0 is the true value 

                contables[r-min_r] << newratio;
                contables[r-min_r] << newratio / true_eigenvalue - 1.;
                contables[r-min_r] << - std::log2( newratio / true_eigenvalue - 1. ); 
                contables[r-min_r] << newratio - true_eigenvalue;
                contables[r-min_r] << std::log2( newratio - true_eigenvalue );
                contables[r-min_r] << u_massnorm;
                contables[r-min_r] << ugrad_massnorm;
                contables[r-min_r] << Float( end - start );
                contables[r-min_r] << nl;

                contables[r-min_r].lg();
            
            }
            
        }

        // estimate the Neumann eigenvalue using one of the recursive estimates 
        if( l == 0 and do_neumann )
        {
        
            LOG << "Estimate PF constant using recursive construction" << nl;
            
            const int num_cells = M.count_tetrahedra();
            const int num_faces = M.count_faces();

            // for later use, collect volumes and diameters of all cells 
            
            std::vector<Float> volumes( num_cells );
            std::vector<Float> diameters( num_cells );

            for( int c = 0; c < num_cells; c++ ) 
            {
                volumes[c]   = M.getMeasure( M.getinnerdimension(), c );
                diameters[c] = M.getDiameter( M.getinnerdimension(), c );
            }

            // for later use, compute the max diameter and the max volume ratio 
            
            Float volumeratio = 1.;
            Float maxdiameter = 0.;
            for( int c1 =    0; c1 < num_cells; c1++ ) 
            for( int c2 = c1+1; c2 < num_cells; c2++ ) 
            {
                volumeratio = maximum( volumeratio, volumes[c1] / volumes[c2] );
                maxdiameter = maximum( maxdiameter, diameters[c1] );
            }
            
            // for later use, compute the xi coefficient
            
            Float Cxi = 0.;
            for( int f = 0; f < num_faces; f++ ) 
            {
                if( M.count_face_tetrahedron_parents(f) != 2 ) continue;
                const auto reflection_matrix_jacobian = M.get_reflection_along_face(f);
                const auto singular_value             = reflection_matrix_jacobian.operator_norm_estimate();
                LOG << singular_value << nl;
                Cxi = maximum( Cxi, singular_value );
            }
            
            //Float Cxi = ( 1. + std::sqrt(6.) * std::sqrt(6.) );
            // Float Cxi = std::sqrt( 1. + power_numerical( std::sqrt(6.) * std::sqrt(6.), 2 ) );
            //Float Cxi = ( 1. + std::sqrt(6.) ) * std::sqrt(6.);

            assert( Cxi <= ( 1. + std::sqrt(6.) ) * std::sqrt(6.) );

            // get the coefficients for the recursion 
            
            Float A  = maxdiameter / std::sqrt(2.);
            Float Ap = maxdiameter / std::sqrt(2.) * std::sqrt( volumeratio ) * Cxi;
            Float B = power_numerical( volumeratio, 1./2. );
            
            LOG << "Constructing the trees. Number of roots: " << num_cells << nl;

            // a tree consists of 1. array of precursors 2. array of rank 3. array of costs 
            // we save the precursors and the costs. The rank will be used in the final computation.

            typedef std::vector<int>   array_prec;
            typedef std::vector<int>   array_rank;
            typedef std::vector<Float> array_cost;
            
            std::vector<array_prec> arrays_of_prec( num_cells, array_prec(num_cells,-1)                                      );
            std::vector<array_rank> arrays_of_rank( num_cells, array_rank(num_cells,-1)                                      );
            std::vector<array_cost> arrays_of_cost( num_cells, array_cost(num_cells, std::numeric_limits<Float>::infinity()) );

            // keep the minimal costs among the trees
            Float minimal_cost_of_tree = std::numeric_limits<Float>::infinity();
            
            // iterate over all the possible roots 
            for( int curr_root = 0; curr_root < num_cells; curr_root++ )
            {
                auto& curr_prec = arrays_of_prec[curr_root];
                auto& curr_rank = arrays_of_rank[curr_root];
                auto& curr_cost = arrays_of_cost[curr_root];
                
                // construct the tree coming from that root 
                {
                    // Priority queue to store {cell, distance} pairs, sorted by distance
                    std::priority_queue<std::pair<int, Float>, std::vector<std::pair<int, Float>>, Compare> pq;

                    // Initialize root data: no precursor, rank zero, some initial costs 
                    curr_prec[curr_root] = -1;
                    curr_rank[curr_root] = 0;
                    curr_cost[curr_root] = A*A;
                    pq.push( { curr_root, curr_cost[curr_root] } );

                    //int security_counter = 0; // if cycles occur    
                    while( not pq.empty() )  // && security_counter++ < num_cells+10
                    {
                        const auto top_entry = pq.top();
                        pq.pop();
                        
                        const int   cell = top_entry.first;
                        const Float cost = top_entry.second;
                        
                        // If this distance is greater than the distance already found, then skip
                        // In this version of Dijkstra, the PQ may have multiple instances of the same node 
                        
                        if( cost > curr_cost[cell] ) continue;

                        // We have retrieved this cell at already the minimum cost and the predecessor is real predecessor
                        // the rank of this cell is the rank of the predecessor + 1
                        
                        if( curr_prec[cell] == -1 ) assert( cell == curr_root );

                        if( curr_prec[cell] != -1 ) curr_rank[cell] = 1 + curr_rank[ curr_prec[cell] ];

                        
                        // Collect all the adjacent cells
                        
                        std::vector<int> adjacent_cells;

                        for( int fi = 0; fi <= M.getinnerdimension(); fi++ )
                        {
                            const int face = M.get_tetrahedron_face( cell, fi );
                            const auto& parents = M.get_tetrahedron_parents_of_face( face );
                            
                            assert( 1 <= parents.size() && parents.size() <= 2 );

                            if( parents.size() == 1 ){
                                assert( parents[0] == cell );
                                continue;
                            }
                            
                            assert( parents[0] == cell or parents[1] == cell );
                            if( parents[0] != cell ) adjacent_cells.push_back( parents[0] );
                            if( parents[1] != cell ) adjacent_cells.push_back( parents[1] );
                        }

                        // Having collected the adjacent cells, iterate over adjacent cells 
                        // compute the costs of each

                        for( int neighbor : adjacent_cells )
                        {
                            // compute the costs of going from the current cell to one its neighbors 
                            
                            auto volume_ratio = maximum( volumes[cell] / volumes[neighbor], volumes[neighbor] / volumes[cell] );
                            
                            Float edge_cost_additive       = volume_ratio; 

                            Float edge_cost_multiplicative = volume_ratio; 
                            
                            // Calculate potential new costs

                            if( curr_rank[cell] == 0 ) {
                                edge_cost_additive       = (Ap+A*B)*(Ap+A*B);
                                edge_cost_multiplicative = 1.; 
                            } else {
                                edge_cost_additive       = (Ap+A*B)*(Ap+A*B) * power_numerical(B*B,curr_rank[cell]);
                                edge_cost_multiplicative = 1.; 
                            }

                            Float new_cost = ( edge_cost_additive + edge_cost_multiplicative * curr_cost[cell] );
                            
                            // If a shorter path is found

                            if( new_cost < curr_cost[neighbor] )
                            {
                                curr_cost[neighbor] = new_cost;
                                curr_prec[neighbor] = cell;  // Set the predecessor
                                pq.push({neighbor, new_cost});
                            }
                        
                        }
                    
                    }
                    //assert( security_counter == num_cells );

                }
                
                // We have constructed the tree!
                
                // DEBUG
                // Check that each cell has been visited (except the root)
                // Either it is the current root or the precursor has is not -1
                for( int i = 0; i < num_cells; i++ ) {
                    assert( i == curr_root || curr_prec[i] != -1 );
                }

                // DEBUG
                // Check that each cell leads to the root and distance is descending along that path to the root
                for( int i = 0; i < num_cells; i++ ) {
                    int c = i;
                    for( int j = 0; j < num_cells && c != curr_root; j++ ) {
                        assert(curr_cost[c] >= curr_cost[curr_prec[c]]); // Ensure non-increasing distances
                        c = curr_prec[c];  // Follow the path to the root
                    }
                    assert(c == curr_root);  // Ensure the cell eventually leads to the root
                }

                // output: total costs of tree
                Float total_costs = std::accumulate( curr_cost.begin(), curr_cost.end(), 0.);
                LOG << "From root " << curr_root << " the total costs are " << total_costs << nl;

                // DEBUG: output the structure of the tree
                // Float total_edge_costs = 0.;
                // for( int i = 0; i < num_cells; i++ ) {
                //     LOG << " At " << i << " prec " << curr_prec[i] << " rank " << curr_rank[i] << " cost " << curr_cost[i] << nl;
                // }

                minimal_cost_of_tree = minimum( minimal_cost_of_tree, total_costs );

            }

            LOG << "minimal costs, std::sqrt : " << minimal_cost_of_tree << space << std::sqrt(minimal_cost_of_tree) << nl;

        }
        
        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

    } 
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
