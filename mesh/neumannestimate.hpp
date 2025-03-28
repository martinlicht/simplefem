#ifndef INCLUDEGUARD_NEUMANNESTIMATE_HPP
#define INCLUDEGUARD_NEUMANNESTIMATE_HPP

#include <algorithm>
#include <limits>
#include <numeric>
#include <queue>
#include <utility>
#include <vector>

#include "../base/include.hpp"
#include "../mesh/mesh.hpp"



// ============================================================================================================
// estimate the Neumann eigenvalue using one of the recursive estimates 
// ============================================================================================================

static Float NeumannEstimate( const Mesh& M ) {


    LOG << "Estimate PF constant using recursive construction" << nl;

    int dim = M.getinnerdimension();
    
    const int num_cells = M.count_simplices(dim);
    const int num_faces = M.count_simplices(dim-1);

    // ======================================================================================
    // for later use, collect volumes and diameters of all cells, as well as reflection norms 
    // ======================================================================================
    
    std::vector<Float> diameters( num_cells, notanumber );
    std::vector<Float> volumes( num_cells, notanumber );
    
    for( int c = 0; c < num_cells; c++ ) 
    {
        diameters[c] = M.getDiameter( dim, c );
        volumes[c]   = M.getMeasure( dim, c );
    }

    std::vector<Float> reflectionnorms( num_faces, notanumber );

    for( int f = 0; f < num_faces; f++ ) 
    {
        if( M.get_supersimplices(dim,dim-1,f).size() != 2 ) continue;
        
        const auto reflection_matrix_jacobian = M.get_reflection_Jacobian_along_face(f);
        const auto singular_value             = reflection_matrix_jacobian.operator_norm_estimate();
        const auto singular_value_inv         = Inverse(reflection_matrix_jacobian).operator_norm_estimate();
        
        reflectionnorms[f] = maximum( singular_value, singular_value_inv );
    }
    
    for( auto v : volumes         ) assert( std::isfinite( v ) );
    for( auto d : diameters       ) assert( std::isfinite( d ) );
    for( auto s : reflectionnorms ) assert( std::isfinite( s ) );
    
    
    // ======================================================================================
    // For each root, a tree consists of 
    // 
    //      1. array of precursors
    //      2. array of rank (level from origin)
    //      3. array of costs 
    // 
    // We save the precursors and the costs. The rank will be used in the final computation.
    // Here, -1 denotes undefined precursor/level, and all costs start with infinity.
    // 
    // For each tree, for each node, we save a coefficient vector.
    // ======================================================================================
    
    LOG << "Constructing the trees. Number of roots: " << num_cells << nl;

    typedef std::vector<int>   array_prec;
    typedef std::vector<int>   array_rank;
    typedef std::vector<Float> array_cost;
    
    std::vector<array_prec> arrays_of_prec( num_cells, array_prec( num_cells, -1                                    ) );
    std::vector<array_rank> arrays_of_rank( num_cells, array_rank( num_cells, -1                                    ) );
    std::vector<array_cost> arrays_of_cost( num_cells, array_cost( num_cells, std::numeric_limits<Float>::infinity()) );

    std::vector<std::vector<FloatVector>> coefficients( num_cells, std::vector<FloatVector>(num_cells, FloatVector(num_cells,0.) ) );

    // Each tree has got a cost, and we keep track fo the minimum cost among the trees found. 

    Float minimal_cost_of_tree = std::numeric_limits<Float>::infinity();
    
    // We iterate over all the possible roots and construct a tree for each.
    for( int curr_root = 0; curr_root < num_cells; curr_root++ )
    {
        // introduce some abbreviations
        auto& curr_prec = arrays_of_prec[curr_root];
        auto& curr_rank = arrays_of_rank[curr_root];
        auto& curr_cost = arrays_of_cost[curr_root];

        auto& curr_coefficients = coefficients[curr_root];
        
        // Comparator for the priority queue (min-heap, based on weights). 
        // Use > to prioritize the pair with the smaller weight
        struct Compare {
            bool operator()(const std::pair<int, Float>& p1, const std::pair<int, Float>& p2) {
                return p1.second > p2.second;
            }
        };

        // Priority queue to store {cell, cost} pairs, sorted by cost
        std::priority_queue<std::pair<int, Float>, std::vector<std::pair<int, Float>>, Compare> pq;

        // ======================================================
        // Initialize root data: 
        // - no precursor,
        // - rank zero, 
        curr_prec[curr_root] = -1;
        curr_rank[curr_root] = 0;

        // - some initial costs, computed like this:
        // We use the Payne-Weinberger bound, or a slight improvement in dimension 2.

        Float Bessel_J11 = 3.83170597020751231561443588630816076656454527428780192876229898991883930951;
        Float natural_poincare_factor = ( dim==2 ? 1./Bessel_J11 : 1./Constants::pi );
        
        curr_coefficients[curr_root][curr_root] = natural_poincare_factor * diameters[curr_root]; // PF constant for convex domains

        curr_cost[curr_root] = curr_coefficients[curr_root].norm_sq();

        // Push the root into the PQ.
        pq.push( { curr_root, curr_cost[curr_root] } );

        // ======================================================
        // Use a Djikstra-type algorithm to compute the queue
        while( not pq.empty() )
        {
            const auto top_entry = pq.top();
            pq.pop();
            
            const int   cell = top_entry.first;
            const Float cost = top_entry.second;
            
            // In this version of Dijkstra, the PQ may have multiple instances of the same node,
            // that is, we queue duplicates rather than update the queue.
            // If this distance is greater than the distance already found, 
            // then this is an outdated entry that we simply skip.
            
            if( cost > curr_cost[cell] ) continue;

            // We have retrieved this cell at already the minimum cost, and the predecessor is real predecessor.
            // The rank of this cell is the rank of the predecessor + 1.
            // The root and only the root has no predecessor.

            LOG << "processing " << cell << " from " << curr_prec[cell] << nl;
            
            if( curr_prec[cell] == -1 ) assert( cell == curr_root );
            if( cell == curr_root ) assert( curr_prec[cell] == -1 );
            
            if( curr_prec[cell] != -1 ) curr_rank[cell] = 1 + curr_rank[ curr_prec[cell] ];

            
            // ============================================================================================================
            // Collect all the adjacent cells and the corresponding connection faces.
            
            std::vector<int> adjacent_cells;
            std::vector<int> connecting_faces;

            for( int fi = 0; fi <= dim; fi++ )
            {
                const int face = M.get_subsimplex( dim, dim-1, cell, fi );
                
                const auto& parents = M.get_supersimplices( dim, dim-1, face );
                
                assert( 1 <= parents.size() && parents.size() <= 2 );

                if( parents.size() == 1 ){
                    assert( parents[0] == cell );
                    continue;
                }
                
                assert( parents[0] == cell or parents[1] == cell );
                assert( parents[0] != cell or parents[1] != cell );
                if( parents[0] != cell ) adjacent_cells.push_back( parents[0] );
                if( parents[1] != cell ) adjacent_cells.push_back( parents[1] );
                
                connecting_faces.push_back(face);
            }

            // LOG << "cell " << cell << " has " << adjacent_cells.size() << " other adjacent cells\n";
            assert( adjacent_cells.size() == connecting_faces.size() );

            
            // ============================================================================================================
            // Having collected the adjacent cells, iterate over adjacent cells. For each, we compute
            // - the cost
            // - the coefficient vector 

            std::vector<Float> adjacent_cells_new_costs( adjacent_cells.size() );

            std::vector<FloatVector> adjacent_cells_new_coeffs( adjacent_cells.size(), FloatVector(num_cells) );

            for( int i = 0; i < adjacent_cells.size(); i++ )
            {
                int neighbor = adjacent_cells[i];

                int face = connecting_faces[i];

                Float volumeratio = volumes[neighbor] / volumes[cell]; // ratio of next cell vs current cell 

                Float reflectionnorm = reflectionnorms[face];

                {
                    // While `reflectionnorm` already contains an estimate, we attempt at improvement
                    // The following is based on estimates in the paper. 

                    assert( cell != neighbor );

                    int fi_index_cell = M.get_subsimplex_index( dim, dim-1, cell,     face );
                    int fi_index_next = M.get_subsimplex_index( dim, dim-1, neighbor, face );

                    int vi_opp_cell = M.get_opposite_subsimplex_index( dim, dim-1, cell,     fi_index_cell );
                    int vi_opp_next = M.get_opposite_subsimplex_index( dim, dim-1, neighbor, fi_index_next );

                    int v_opp_cell = M.get_subsimplex( dim, 0, cell,     vi_opp_cell );
                    int v_opp_next = M.get_subsimplex( dim, 0, neighbor, vi_opp_next );
                    
                    auto h_cell = M.getHeightVector( dim, cell,     vi_opp_cell );
                    auto h_next = M.getHeightVector( dim, neighbor, vi_opp_next );
                    
                    Float a = h_cell.norm() / h_next.norm();
                    
                    Assert( is_numerically_close( a, 1./volumeratio ), a, 1./volumeratio );
                
                    auto z_cell = M.getCoordinates().getdata_by_vertex( v_opp_cell );
                    auto z_next = M.getCoordinates().getdata_by_vertex( v_opp_next );

                    auto b_cell = z_cell - ( z_cell * h_cell ) / h_cell.norm_sq() * h_cell;
                    auto b_next = z_next - ( z_next * h_next ) / h_next.norm_sq() * h_next;

                    Float c = ( b_cell - b_next ).norm() / h_next.norm();

                    Float facebased_singular_max = 0.5 * std::sqrt( square( a + 1. ) + c ) + 0.5 * std::sqrt( square( a - 1. ) + c );

                    LOG << "Reflection estimates " << reflectionnorm << space << facebased_singular_max << nl;

                    reflectionnorm = minimum( reflectionnorm, facebased_singular_max );

                }

                // We estimate the Poincare constant with boundary conditions on the simplex. 
                // With that, we build the coefficients in the recursion. 

                Float pf_simplex_with_bc = minimum( natural_poincare_factor, 1./sqrt(2.) ) * diameters[cell];
                
                Float  A_via_bc = pf_simplex_with_bc;
                Float Ap_via_bc = sqrt(volumeratio) * reflectionnorm * pf_simplex_with_bc;

                Float A =  A_via_bc; // coefficient of new cell 

                Float B = Ap_via_bc; // coefficient of previous cell 

                Float C = sqrt( volumeratio ); // this is the recursive coefficient

                // LOGPRINTF( "SCALAR nat=%f diam=%f A=%f B=%f C=%f \n", natural_poincare_factor, diameters[cell],  A, B, C );

                FloatVector new_coeffs = C * curr_coefficients[cell];

                new_coeffs[cell] += B;

                new_coeffs[neighbor] += A;

                Float new_cost = new_coeffs.norm_sq();

                // Finally, we set the costs and coefficient vectors for this particular cell.

                adjacent_cells_new_costs[i] = new_cost;

                adjacent_cells_new_coeffs[i] = new_coeffs;
            }

            for( auto t : adjacent_cells_new_coeffs ) assert( t.is_finite() );

            // Having found the costs, do the relevant updates 
                
            for( int i = 0; i < adjacent_cells.size(); i++ )
            { 
                int neighbor = adjacent_cells[i];

                Float new_cost = adjacent_cells_new_costs[i];

                auto new_coeffs = adjacent_cells_new_coeffs[i];

                // If a shorter path is found

                if(- new_cost < curr_cost[neighbor] )
                {
                    // LOG << "update: " << curr_cost[neighbor] << " to " << new_cost << nl;
                    
                    curr_cost[neighbor] = new_cost;
                
                    curr_prec[neighbor] = cell;  

                    curr_coefficients[neighbor] = new_coeffs;

                    pq.push({neighbor, new_cost});
                }
            }
        }
        
        // =========================================================================================
        // We have constructed the complete tree!
        
        // DEBUG
        // Check that each cell has been visited (except the root)
        // Either it is the current root or the precursor has is not -1
        // Additionally, check that the rank makes sense
        for( int i = 0; i < num_cells; i++ ) {
            
            if( i == curr_root ) 
            {
                assert( curr_prec[i] == -1 );
                assert( curr_rank[i] == 0 );
                continue;
            };
            
            assert( curr_prec[i] != -1 );
            assert( 0 <= curr_prec[i] && curr_prec[i] < num_cells );
            assert( curr_rank[i] == 1 + curr_rank[curr_prec[i]] );
        }

        // DEBUG
        // Check that each cell has costs corresponding to the coefficient vectors l2 norm squared
        for( int i = 0; i < num_cells; i++ ) {
            Float norm = curr_coefficients[i].norm_sq();
            Assert( curr_cost[i] == norm, curr_cost[i], norm );
        }

        // DEBUG
        // Check that each cell leads to the root and distance is descending along that path to the root
        for( int i = 0; i < num_cells; i++ ) {
            int c = i;
            for( int j = 0; j < num_cells && c != curr_root; j++ ) {
                // Assert( curr_cost[c] >= curr_cost[curr_prec[c]], curr_cost[c], curr_cost[curr_prec[c]], c, curr_prec[c] ); // Ensure non-increasing distances
                c = curr_prec[c];  // Follow the path to the root
            }
            assert(c == curr_root);  // Ensure the cell eventually leads to the root
        }

        // ================================================================================================================
        // We compute the total costs of the tree. If necessary, we update the minimal cost among all trees found.
        Float total_costs = std::accumulate( curr_cost.begin(), curr_cost.end(), 0.);
        
        LOG << "From root " << curr_root << " the total costs are " << total_costs << " -> " << sqrt(total_costs) << nl;

        // DEBUG: output the structure of the tree
        for( int i = 0; i < num_cells; i++ ) {
            LOG << " At " << i << " prec " << curr_prec[i] << " rank " << curr_rank[i] << " cost " << curr_cost[i] << nl;
        }
        for( int i = 0; i < num_cells; i++ ) {
            LOG << curr_coefficients[i] << nl;
        }
        
        minimal_cost_of_tree = minimum( minimal_cost_of_tree, total_costs );

    }

    LOG << "minimal costs, sqrt : " << minimal_cost_of_tree << space << std::sqrt(minimal_cost_of_tree) << nl;

    return minimal_cost_of_tree;
}




#endif
