#ifndef INCLUDEGUARD_SPANNING_TREE_LISTER_HPP
#define INCLUDEGUARD_SPANNING_TREE_LISTER_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <utility> // For std::pair

#include "../base/include.hpp"
#include "mesh.hpp"






// ============================================================================================================
// Auxiliary function to generate all strictly ascending combinations of k elements from 0 to N-1.
// ============================================================================================================

std::vector<std::vector<int>> generate_combinations( int N, int k );

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

    assert( result.size() == binomial_integer_secured( N, k ) );

    return result;
}






// ========================================================================================================================
// The following is numerous auxiliary functions that address the interleaving of two lists. 
// 
// In other words, given two totally ordered sets, list all possible total orders on the union of those two sets
// that include the original two total orders. 
// ========================================================================================================================

std::vector<std::vector<int>> interleaveLists( 
    const std::vector<int>& list1,
    const std::vector<int>& list2
);

void interleaveLists(
    const std::vector<int>& list1,
    const std::vector<int>& list2,
    std::vector<int>& current,
    int index1,
    int index2,
    std::vector<std::vector<int>>& results
);

std::vector<std::vector<int>> interleaveLists( 
    const std::vector<std::vector<int>>& listlists1,
    const std::vector<std::vector<int>>& listlists2
);

std::vector<std::vector<int>> interleaveLists( 
    const std::vector<std::vector<std::vector<int>>>& listlistlists
);

// ============================================================================================================
// Given two lists of integers, list1 and list2, create all possible interleavings of those lists.
// ============================================================================================================
std::vector<std::vector<int>> interleaveLists( 
    const std::vector<int>& list1,
    const std::vector<int>& list2
){
    std::vector<std::vector<int>> results;
    std::vector<int> current;
    interleaveLists( list1, list2, current, 0, 0, results );
    return results;
}

// ============================================================================================================
// Internal auxiliary function function to interleave list1 and list2 using backtracking
// ============================================================================================================
void interleaveLists(
    const std::vector<int>& list1,
    const std::vector<int>& list2,
    std::vector<int>& current,
    int index1,
    int index2,
    std::vector<std::vector<int>>& results
){

    // If we've used all elements from both lists, store the result
    if (index1 == list1.size() && index2 == list2.size()) {
        results.push_back(current);
        return;
    }

    // If there are remaining elements in list1, add the next one
    if (index1 < list1.size()) {
        current.push_back(list1[index1]);
        interleaveLists(list1, list2, current, index1 + 1, index2, results);
        current.pop_back(); // Backtrack
    }

    // If there are remaining elements in list2, add the next one
    if (index2 < list2.size()) {
        current.push_back(list2[index2]);
        interleaveLists(list1, list2, current, index1, index2 + 1, results);
        current.pop_back(); // Backtrack
    }
}

// ============================================================================================================
// Given two lists of lists of integers, listlist1 and listlist2, create all possible interleavings 
// of pairs list1 and list2 from those lists of lists.
// ============================================================================================================
std::vector<std::vector<int>> interleaveLists( 
    const std::vector<std::vector<int>>& listlists1,
    const std::vector<std::vector<int>>& listlists2
){
    std::vector<std::vector<int>> results;

    for( const auto& list1 : listlists1 )
    for( const auto& list2 : listlists2 )
    {
        const auto temp = interleaveLists( list1, list2 );
        results.insert( results.end(), temp.begin(), temp.end() );
    }

    return results;
}

// ============================================================================================================
// Given a list of S list of lists of integers, listlistlists, 
// list all the possible interleavings of the lists of lits of integers 
// listlistlists[0], listlistlists[1], ..., listlistlists[S].
// ============================================================================================================
std::vector<std::vector<int>> interleaveLists( 
    const std::vector<std::vector<std::vector<int>>>& listlistlists
){
    if( listlistlists.size() == 0 ) return {};

    std::vector<std::vector<int>> result = listlistlists[0];
    
    for( int t = 1; t < listlistlists.size(); t++ )
    {
        result = interleaveLists( result, listlistlists[t] );
    }

    // LOG << "LLL+" << result.size() << nl;
    
    return result;
}





// ==========================================================================================================================
// Functions to create spanning trees.
// ==========================================================================================================================

struct Link{
    int index;
    int first;
    int second;
};

std::vector<Link> get_nodes_of_links( const Mesh& mesh );

bool check_spanning_tree( 
    const std::vector<Link>& nodes_of_links, 
    const std::vector<int>& tree_candidate 
);

std::pair< std::vector<int>, std::vector<std::vector<int>> > list_face_spanning_trees( const Mesh& mesh );


// ==========================================================================================================================
// Given a mesh, return the list of face-links. This provides the basic structure to discuss face-spanning sets.
// ==========================================================================================================================
std::vector<Link> get_nodes_of_links( const Mesh& mesh )
{
    // check that input mesh is reasonable and gather basic data 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const int num_volumes = mesh.count_simplices(dim);
    const int num_faces   = mesh.count_simplices(dim-1);

    assert( num_faces >= num_volumes-1 );

    // enumerate all the links' nodes

    std::vector<Link> nodes_of_links;
    nodes_of_links.reserve( num_faces );
    
    for( int f = 0; f < num_faces; f++ ) 
    {
        const auto parent_volumes = mesh.get_supersimplices( dim, dim-1, f );
        
        if( parent_volumes.size() == 1 ) continue; 
        
        Link link = { .index = f, .first = parent_volumes[0], .second = parent_volumes[1] };
        assert( 0 <= link.first  && link.first  < num_volumes );
        assert( 0 <= link.second && link.second < num_volumes );
        assert( link.first != link.second );
        
        nodes_of_links.push_back( link );
    }
    
    return nodes_of_links;
}

// ==========================================================================================================================
// Check whether a given subgraph is a spanning tree.
// ==========================================================================================================================
bool check_spanning_tree( 
        const std::vector<Link>& nodes_of_links, 
        const std::vector<int>& tree_candidate 
)
{
    for( const auto& link : nodes_of_links ) assert( link.first >= 0 && link.second >= 0 );
    for( const auto& link : nodes_of_links ) assert( link.index >= 0 );

    for( const auto& link : tree_candidate ) assert( 0 <= link && link < nodes_of_links.size() );

    int number_of_nodes = 0; 
    for( const auto& link : nodes_of_links ) 
    {
        number_of_nodes = std::max( number_of_nodes, 1 + std::max( link.first, link.second ) );
    }

    if( tree_candidate.size() != number_of_nodes - 1 ) return false;

    // check whether all indices between 0 and number_of_nodes are present in one pair of the tree candidate, 
    // which means the tree candidate covers all nodes 
    std::vector<int> visited( number_of_nodes, false );
    for( const auto& l : tree_candidate ) 
    {
        const auto& link = nodes_of_links[l];
        int n1 = link.first;
        int n2 = link.second;
        assert( 0 <= n1 && n1 < number_of_nodes );
        assert( 0 <= n2 && n2 < number_of_nodes );
        visited[n1] = true;
        visited[n2] = true;
    }
    for( const bool v : visited ) if( not v ) return false;

    // create the adjaceny list for all nodges
    std::vector<std::vector<int>> adjacency_list(number_of_nodes);
    for( const auto& l : tree_candidate) {
        const auto& link = nodes_of_links[l];
        int u = link.first;
        int v = link.second;
        adjacency_list[u].push_back(v);
        adjacency_list[v].push_back(u);
    }

    // traverse the tree candidate via DFS and check whether all links have been visited
    std::fill( visited.begin(), visited.end(), false );
    std::vector<int> S;
    S.reserve( adjacency_list.size() );
    S.push_back(0); // Start DFS from node 0
    visited[0] = true;

    int visited_count = 0;
    while( not S.empty() ) 
    {
        int node = S.front();
        S.pop_back();
        visited_count++;

        for( int neighbor : adjacency_list[node] )
        {
            if( not visited[neighbor] )
            {
                visited[neighbor] = true;
                S.push_back(neighbor);
            }
        }
    }

    // The graph is connected if all nodes have been visited
    return visited_count == number_of_nodes;
}

// ==========================================================================================================================
// List all the face spanning trees.
// ==========================================================================================================================
std::pair< std::vector<int>, std::vector<std::vector<int>> > list_face_spanning_trees( const Mesh& mesh )
{
    // Check that input mesh is reasonable.
    
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const int num_volumes = mesh.count_simplices(dim);
    const int num_faces   = mesh.count_simplices(dim-1);

    assert( num_faces >= num_volumes-1 );

    // Enumerate all the links' nodes and generate the possible candidates for spanning trees.
    // The construction is brute force: we list all subgraphs whose edge-cardinality satisfies 
    // the necessary criterion for being an spanning tree.
    
    const std::vector<Link> nodes_of_links = get_nodes_of_links( mesh );
    
    LOG << "generate combinations: " << nodes_of_links.size() << space << num_volumes << nl;
    
    auto spanning_tree_candidates = generate_combinations( nodes_of_links.size(), num_volumes-1 );

    // LOG << nodes_of_links.size() << space << binomial_integer( nodes_of_links.size(), num_volumes-1 ) << nl;
    // for( const auto& candidate : spanning_tree_candidates ) {
    //     for( const int link : candidate ) LOG << e << space;
    //     LOG << "---" << nl;
    // }

    LOG << "check for trees: " << spanning_tree_candidates.size() << nl;

    std::vector<bool> is_tree( spanning_tree_candidates.size() );
    
    for( int t = 0; t < spanning_tree_candidates.size(); t++ )
    {
        const auto& candidate = spanning_tree_candidates[t];
        
        is_tree[t] = check_spanning_tree( nodes_of_links, candidate );
    }

    std::vector<std::vector<int>> spanning_trees;
    spanning_trees.reserve( spanning_tree_candidates.size() );

    for( int t = 0; t < spanning_tree_candidates.size(); t++ )
    {
        if( is_tree[t] ) 
            spanning_trees.push_back( spanning_tree_candidates[t] );
    }
    spanning_trees.shrink_to_fit();

    std::vector<int> index2face( nodes_of_links.size() );
    for( int i = 0; i < nodes_of_links.size(); i++ ) index2face[i] = nodes_of_links[i].index;

    return { index2face, spanning_trees };
}




















// ============================================================================================================
// Functionality to list the topological sortings of a given spanning tree. 
// ============================================================================================================

std::vector<std::vector<int>> list_ordered_face_spanning_trees( 
    const Mesh& mesh, 
    // const std::vector<int>& index2face, 
    const std::vector<int>& spanning_tree
);

std::vector<std::vector<int>> list_topological_node_sortings_of_subtree(
    const Mesh& mesh, 
    const std::vector<std::vector<int>>& adjacency_list, 
    int forbidden_node,
    int current_node
);
    
// ============================================================================================================
// List the topological sortings of a given spanning tree. 
// Here, the mesh and adjacency list (neighbors of each node) are already known. We start at the `current_node`
// and build from there. The `forbidden_node` effectively marks the predecessor and blocks some parts of the 
// tree from re-exploration.
// ============================================================================================================
std::vector<std::vector<int>> list_topological_node_sortings_of_subtree(
    const Mesh& mesh, 
    const std::vector<std::vector<int>>& adjacency_list, 
    int forbidden_node,
    int current_node
){
    // check consitency and gather data 

    assert( forbidden_node == Mesh::nullindex || ( 0 <= forbidden_node && forbidden_node < adjacency_list.size() ) );
    assert( 0 <= current_node && current_node < adjacency_list.size() );
    assert( adjacency_list.size() == mesh.count_simplices( mesh.getinnerdimension() ) );
    for( const auto& vs : adjacency_list ) for( const auto& v : vs ) assert( 0 <= v && v < mesh.count_simplices( mesh.getinnerdimension() ) );
    
    const std::vector<int>& start_neighbors = adjacency_list[current_node];

    if( forbidden_node != Mesh::nullindex ) 
        assert( std::find( start_neighbors.cbegin(), start_neighbors.cend(), forbidden_node ) != start_neighbors.end() );
    
    // special cases that signal the end of the recursion
    
    const int number_of_neighbors = start_neighbors.size();

    if( number_of_neighbors == 0 ){
        assert( forbidden_node == Mesh::nullindex );
        return {{current_node}};
    }

    if( number_of_neighbors == 1 && forbidden_node != Mesh::nullindex ){
        return {{current_node}};
    }

    // Recursion: list the topological sortings coming the successors

    std::vector<std::vector<std::vector<int>>> top_sorts_of_succs;
    top_sorts_of_succs.reserve( number_of_neighbors );

    for( int t = 0; t < number_of_neighbors; t++ )
    {
        int succ = start_neighbors[t];

        if( succ == forbidden_node ) continue;

        const auto& current = list_topological_node_sortings_of_subtree( mesh, adjacency_list, current_node, succ );

        assert( current.size() > 0 );

        top_sorts_of_succs.push_back( current );
    }

    assert( top_sorts_of_succs.size() > 0 );

    // Interleave the successors and prepend the current source 

    std::vector<std::vector<int>> results;
    
    results = interleaveLists( top_sorts_of_succs );

    for( auto& list : results ) list.insert( list.begin(), current_node );
    
    assert( results.size() > 0 );
    
    return results;

}

// =============================================================================
// Given a face-graph spanning tree for a mesh, 
// together with an index of the faces, 
// produce all topological sortings of that spanning tree. 
// =============================================================================
std::vector<std::vector<int>> list_ordered_face_spanning_trees( 
    const Mesh& mesh, 
    // const std::vector<int>& index2face, 
    const std::vector<int>& spanning_tree
){
    // Check that input mesh is reasonable and collect basic data 

    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const int num_volumes = mesh.count_simplices(dim);
    const int num_faces   = mesh.count_simplices(dim-1);

    assert( num_faces >= num_volumes-1 );

    const int number_of_nodes = num_volumes;

    // Enumerate all the links' nodes
    
    const std::vector<Link> nodes_of_links = get_nodes_of_links( mesh );
    
    std::vector<std::vector<int>> results;
    
    // Create the node adjacency list for the spanning tree 

    std::vector<std::vector<int>> adjacency_list(number_of_nodes);
    for( const auto& l : spanning_tree ) 
    {
        int u = nodes_of_links[l].first;
        int v = nodes_of_links[l].second;
        adjacency_list[u].push_back(v);
        adjacency_list[v].push_back(u);
    }

    // create the link adjacency list for the tree subgraph 
    // std::vector<std::vector<int>> adjacency_list(num_faces);
    // for( const auto& e1 : spanning_tree ) 
    // for( const auto& e2 : spanning_tree ) 
    // {
    //     int u1 = nodes_of_links[e1].first;
    //     int v1 = nodes_of_links[e1].second;
    //     int u2 = nodes_of_links[e2].first;
    //     int v2 = nodes_of_links[e2].second;
    //     if( u1 == u2 and v1 == v2 ) continue;
    //     if( u1 != u2 and v1 != v2 ) continue;
    //     adjacency_list[e1].push_back(e2);
    //     adjacency_list[e2].push_back(e1);
    // }

    // For each node, create all rooted enumerations of the tree
    
    for( int v = 0; v < number_of_nodes; v++ )
    {
        const auto& list = list_topological_node_sortings_of_subtree( mesh, adjacency_list, Mesh::nullindex, v );

        results.insert( results.end(), list.begin(), list.end() );
    }
        
    return results;
}



#endif // INCLUDEGUARD_SPANNING_TREE_LISTER_HPP
