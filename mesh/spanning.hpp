#ifndef SPANNING_TREE_LISTER_HPP
#define SPANNING_TREE_LISTER_HPP

#include <vector>
#include <queue>
#include <utility> // For std::pair

#include "../basic.hpp"
#include "mesh.hpp"

////////////////////////////////////////
////////  forward declarations /////////
////////////////////////////////////////

struct Edge{
    int index;
    int first;
    int second;
};

std::vector<std::vector<int>> generate_combinations( int N, int k );

bool check_spanning_tree( 
        const std::vector<Edge>& nodes_of_edges, 
        const std::vector<int>& tree_candidate 
);

std::pair< std::vector<int>, std::vector<std::vector<int>> > list_face_spanning_trees( const Mesh& mesh );

////////////////////////////////////////





// Function to generate all combinations of k elements from 0 to N-1
std::vector<std::vector<int>> generate_combinations( int N, int k ){
    
    std::vector<std::vector<int>> result;
    std::vector<int> combination(k);
    
    assert( 0 <= k and k < N );
    
    // Initialize the combination vector with the first k indices
    for( int i = 0; i < k; ++i ) 
    {
        combination[i] = i;
    }

    // Generate combinations
    while( true )
    {
        
        result.push_back(combination);

        // Find the first element from the right that can be incremented
        int i = k - 1;
        while( i >= 0 && combination[i] == N - k + i )
        {
            --i;
        }

        // If no such element exists, we are done
        if( i < 0 )
        {
            break;
        }

        // Increment this element
        ++combination[i];

        // Adjust all subsequent elements
        for( int j = i + 1; j < k; ++j )
        {
            combination[j] = combination[j - 1] + 1;
        }
    }

    assert( result.size() == binomial_integer( N, k ) );

    return result;
}






bool check_spanning_tree( 
        const std::vector<Edge>& nodes_of_edges, 
        const std::vector<int>& tree_candidate 
)
{
    for( const auto& e : nodes_of_edges ) assert( e.first >= 0 and e.second >= 0 );
    for( const auto& e : nodes_of_edges ) assert( e.index >= 0 );

    for( const auto& e : tree_candidate ) assert( 0 <= e and e < nodes_of_edges.size() );

    int number_of_nodes = 0; 
    for( const auto& e : nodes_of_edges ) 
    {
        number_of_nodes = std::max( number_of_nodes, 1 + std::max( e.first, e.second ) );
    }

    if( tree_candidate.size() != number_of_nodes - 1 ) return false;

    // check whether all indices between 0 and number_of_nodes are present in one pair of the tree candidate, 
    // which means the tree candidate covers all nodes 
    std::vector<int> visited( number_of_nodes, false );
    for( const auto& e : tree_candidate ) 
    {
        int n1 = nodes_of_edges[e].first;
        int n2 = nodes_of_edges[e].second;
        assert( 0 <= n1 and n1 < number_of_nodes );
        assert( 0 <= n2 and n2 < number_of_nodes );
        visited[n1] = true;
        visited[n2] = true;
    }
    for( const bool v : visited ) if( not v ) return false;



    // create the adjaceny list for all nodges
    std::vector<std::vector<int>> adjacency_list(number_of_nodes);
    for (const auto& e : tree_candidate) {
        int u = nodes_of_edges[e].first;
        int v = nodes_of_edges[e].second;
        adjacency_list[u].push_back(v);
        adjacency_list[v].push_back(u);
    }

    


    // traverse the tree candidate via BFS and check whether all edges have been visited
    std::fill( visited.begin(), visited.end(), false );
    std::queue<int> q;
    q.push(0); // Start BFS from node 0
    visited[0] = true;

    int visited_count = 0;
    while( not q.empty() ) 
    {
        int node = q.front();
        q.pop();
        visited_count++;

        for( int neighbor : adjacency_list[node] )
        {
            if( not visited[neighbor] )
            {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }

    // The graph is connected if all nodes have been visited
    return visited_count == number_of_nodes;
    
}





std::pair< std::vector<int>, std::vector<std::vector<int>> > list_face_spanning_trees( const Mesh& mesh )
{
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const int num_volumes = mesh.count_simplices(dim);
    const int num_faces   = mesh.count_simplices(dim-1);

    assert( num_faces >= num_volumes-1 );

    // generate the list of parents of each face with two parents 
    std::vector<Edge> nodes_of_edges;
    nodes_of_edges.reserve( num_faces );
    
    for( int f = 0; f < num_faces; f++ ) 
    {
        const auto parent_volumes = mesh.getsupersimplices( dim, dim-1, f );
        
        if( parent_volumes.size() == 1 ) continue; 
        
        Edge edge = { .index = f, .first = parent_volumes[0], .second = parent_volumes[1] };
        assert( 0 <= edge.first  and edge.first  < num_volumes );
        assert( 0 <= edge.second and edge.second < num_volumes );
        assert( edge.first != edge.second );
        
        nodes_of_edges.push_back( edge );
    }
    
    // generate the ordered spanning tree candidates and filter them out 
    auto spanning_tree_candidates = generate_combinations( nodes_of_edges.size(), num_volumes-1 );

    // LOG << nodes_of_edges.size() << space << binomial_integer( nodes_of_edges.size(), num_volumes-1 ) << nl;
    // for( const auto& candidate : spanning_tree_candidates ) {
    //     for( const int e : candidate ) LOG << e << space;
    //     LOG << "---" << nl;
    // }

    for( int t = 0; t < spanning_tree_candidates.size(); )
    {
        const auto& candidate = spanning_tree_candidates[t];

        bool is_tree = check_spanning_tree( nodes_of_edges, candidate );

        if( not is_tree ) {
            spanning_tree_candidates.erase( spanning_tree_candidates.begin() + t );
        } else {
            t++;
        }
    }

    std::vector<int> index2face( nodes_of_edges.size() );
    for( int i = 0; i < nodes_of_edges.size(); i++ ) index2face[i] = nodes_of_edges[i].index;

    return { index2face, spanning_tree_candidates };

}

#endif // SPANNING_TREE_LISTER_HPP
