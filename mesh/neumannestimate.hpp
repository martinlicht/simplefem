#ifndef INCLUDEGUARD_NEUMANNESTIMATE_HPP
#define INCLUDEGUARD_NEUMANNESTIMATE_HPP
















#include <algorithm>
#include <numeric>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "../../basic.hpp"
#include "../../mesh/mesh.hpp"



// TODO: can this be a local lambda?
// TODO: Double check correctness
// Comparator for the priority queue (min-heap, based on weights). Use > to prioritize the pair with the smaller weight
struct Compare {
    bool operator()(const std::pair<int, Float>& p1, const std::pair<int, Float>& p2) {
        return p1.second > p2.second;
    }
};



Float NeumannEstimate( const Mesh& M ) {

    // estimate the Neumann eigenvalue using one of the recursive estimates 

    LOG << "Estimate PF constant using recursive construction" << nl;

    int dim = M.getinnerdimension();
    
    const int num_cells = M.count_simplices(dim);
    const int num_faces = M.count_simplices(dim-1);

    // for later use, collect volumes and diameters of all cells 
    
    std::vector<Float> volumes( num_cells );
    std::vector<Float> diameters( num_cells );

    for( int c = 0; c < num_cells; c++ ) 
    {
        volumes[c]   = M.getMeasure( dim, c );
        diameters[c] = M.getDiameter( dim, c );
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
        if( M.get_supersimplices(dim,dim-1,f).size() != 2 ) continue;
        const auto reflection_matrix_jacobian = IdentityMatrix(1); // TODO: Replace this computation with something more universal
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
    Float minimal_cost_of_tree = numeric_limits<Float>::infinity();
    
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




#endif