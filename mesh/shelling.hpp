#ifndef SHELLING_LISTER_HPP
#define SHELLING_LISTER_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <utility> // For std::pair

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "mesh.hpp"




std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh
);

std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list
);

void generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list,
    std::vector<std::vector<int>>& shellings_found,
    const std::vector<int> current_prefix,
    std::vector<int> remaining_nodes
);




std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh
){
    const int dim = mesh.getinnerdimension();

    std::vector<bool> acceptable_face_list( mesh.count_simplices(dim-1) );
    
    for( int face = 0; face < mesh.count_simplices(dim-1); face++ ) 
    {
        const auto& parents = mesh.getsupersimplices(dim,dim-1,face);
        acceptable_face_list[face] = ( parents.size() != 1 );
    }

    return generateShellings( mesh, acceptable_face_list );
}


std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list
){
    const int dim = mesh.getinnerdimension();

    assert( acceptable_face_list.size() == mesh.count_simplices( dim-1 ) );

    std::vector<std::vector<int>> shellings_found;
    
    std::vector<int> current_prefix;

    std::vector<int> remaining_nodes( mesh.count_simplices( dim ) );
    for( int node = 0; node < remaining_nodes.size(); node++ ) remaining_nodes[node] = node;

    generateShellings( mesh, acceptable_face_list, shellings_found, current_prefix, remaining_nodes );
    
    return shellings_found;
}




void generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list,
    std::vector<std::vector<int>>& shellings_found,
    const std::vector<int> current_prefix,
    std::vector<int> remaining_nodes
){
    
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );


    // count all the simplices of the mesh 
    const auto counts = mesh.count_simplices();


    // if enough shellings have been found already, then return 
    if( false and shellings_found.size() >= 1 )
        return;


    // if the current prefix is already containing all volumes, 
    // then a shelling has been found and can be added 
    if( current_prefix.size() == counts[dim] )
    {
        shellings_found.push_back( current_prefix );
        LOG << "found" << nl;
        return;
    }


    // There are still volumes to be added
    // Run over all volumes and check whether they can be added 
    for( int node : remaining_nodes )
    {
        
        assert( not contains( current_prefix, node ) );
        
        
        {
            for( int t : current_prefix ) LOG << t << space;
            LOG << node << nl;    
        }
        

        // list the faces adjacent to the node
        // and check whether it is connected to one of th prefix nodes  
        
        const auto faces_of_node = mesh.getsubsimplices( dim, dim-1, node ).getvalues();
        assert( faces_of_node.size() == dim+1 );
        std::vector<bool> face_is_connected( dim+1, false );

        for( int index_f = 0; index_f < faces_of_node.size(); index_f++ )
        {
            const int face = faces_of_node[index_f];
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( 1 <= parents.size() and parents.size() <= 2 );

            if( parents.size() == 1 ) continue;

            assert( parents[0] == node or parents[1] == node );
            assert( parents[0] != node or parents[1] != node );

            if( parents[0] != node and contains( current_prefix, parents[0] ) ) face_is_connected[index_f] = true;
            if( parents[1] != node and contains( current_prefix, parents[1] ) ) face_is_connected[index_f] = true;
        }

        
        // OPTIMIZATION
        // if current_prefix is not empty and volume node has no common face with the volumes in current_prefix, then skip 
        
        if( current_prefix.size() != 0 )
        if( not std::any_of( face_is_connected.begin(), face_is_connected.end(), [](bool b) { return b; } ) )
            continue;


        // If the connection to the previous volumes is not via an acceptable face, then skip
        
        bool acceptable_connection = false;
        
        for( int index_f = 0; index_f <= dim; index_f++ ){
            
            const auto face = faces_of_node[index_f];
            
            if( acceptable_face_list[face] and ( current_prefix.size() == 0 or face_is_connected[index_f] ) )
                acceptable_connection = true;
        }

        if( not acceptable_connection ) {
            LOG << "no acceptable face" << nl;
            continue;
        }
        
        
        // for each proper subsimplex index_f of the volume s
        // if index_f is a subsimplex of one of the previous volumes 
        // then it must be contained in one of the connected faces 

        bool is_compatible = true;

        // List all proper subsimplices of the node, excluding itself and its faces 
        for( int d = 0; d <= dim-2; d++ )
        {
            if( not is_compatible ) break;
            
            const auto subsimplices_of_node = mesh.getsubsimplices( dim, d, node ).getvalues();

            // For each such proper subsimplex of the node ...
            for( const auto sub : subsimplices_of_node )
            {
                
                if( not is_compatible ) break;

                // check whether it has got a supersimplex within the prefix
                const auto supersimplices_of_sub = mesh.getsupersimplices( dim, d, sub );

                bool sub_of_prefix = std::any_of( 
                    supersimplices_of_sub.begin(), supersimplices_of_sub.end(),
                    [&]( int val ) { return contains<int>( current_prefix, val ); } 
                );

                // if not, there is nothing to check 
                if( not sub_of_prefix ) continue;

                // otherwise, check whether its contained within one of the active faces 
                bool within_connected_face = false;

                const auto faces_containing_sub = mesh.getsupersimplices( dim-1, d, sub );

                for( int index_f = 0; index_f < faces_of_node.size(); index_f++ ){
                    if( not face_is_connected[index_f] ) continue;
                    if( contains( faces_containing_sub, faces_of_node[index_f] ) )
                        within_connected_face = true;
                }

                // if sub is not within a connected face, then node is not compatible 
                if( not within_connected_face )
                    is_compatible = false;

            }
        }

        // We have determined whether node is compatible. If not, then continue with the next one

        // is_compatible = true; // TODO

        if( not is_compatible ) {
            LOG << "reject" << nl;
            continue;
        }

        // node is compatible with the prefix, then recursion 

        {
            auto next_prefix = current_prefix;
            next_prefix.reserve( counts[dim] );
            next_prefix.push_back( node );

            auto next_remaining_nodes = remaining_nodes;
            auto it = std::find( next_remaining_nodes.begin(), next_remaining_nodes.end(), node );
            assert( it != next_remaining_nodes.end() );
            next_remaining_nodes.erase( it );

            generateShellings( mesh, acceptable_face_list, shellings_found, next_prefix, next_remaining_nodes );
        }

    }
   
}























std::vector<std::vector<std::pair<int,Float>>> generate_ranked_shelling(
    const Mesh& mesh
);

void generate_ranked_shelling(
    const Mesh& mesh,
    std::vector<std::vector<std::pair<int,Float>>>& shellings_found,
    std::vector<std::pair<int,Float>> current_prefix,
    std::vector<std::pair<int,Float>> remaining_nodes
);


std::vector<std::vector<std::pair<int,Float>>> generate_ranked_shelling(
    const Mesh& mesh
){
    const int dim = mesh.getinnerdimension();

    std::vector<std::vector<std::pair<int,Float>>> ret;

    for( int t = 0; t < mesh.count_simplices(dim); t++ ) 
    {
        std::vector<std::vector<std::pair<int,Float>>> shellings_found;
        
        std::vector<std::pair<int,Float>> current_prefix;

        std::vector<std::pair<int,Float>> remaining_nodes;
        remaining_nodes.reserve( mesh.count_simplices( dim ) );
        for( int s = 0; s < mesh.count_simplices(dim); s++ ) 
            if( s != t ) 
                remaining_nodes.push_back( { s, std::numeric_limits<Float>::infinity() } );

        current_prefix.push_back( { t, 0. } );

        // CHECK that current and remaining nodes are disjoint
        for( const auto& node1 : current_prefix ) 
        for( const auto& node2 : remaining_nodes ) 
            assert( node1.first != node2.first );


        generate_ranked_shelling( mesh, shellings_found, current_prefix, remaining_nodes );
        
        ret.insert( ret.end(), shellings_found.begin(), shellings_found.end() ); 
    }
    return ret;
}


void generate_ranked_shelling(
    const Mesh& mesh,
    std::vector<std::vector<std::pair<int,Float>>>& shellings_found,
    std::vector<std::pair<int,Float>> current_prefix,
    std::vector<std::pair<int,Float>> remaining_nodes
){
    
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );


    // CHECK that current and remaining nodes are disjoint
    for( const auto& node1 : current_prefix ) 
    for( const auto& node2 : remaining_nodes ) 
        assert( node1.first != node2.first );

    
    // count all the simplices of the mesh 
    const auto counts = mesh.count_simplices();

    
    // if enough shellings have been found already, then return 
    if( false and shellings_found.size() >= 1 )
        return;


    // if the current prefix is already containing all volumes, 
    // then a shelling has been found and can be added 
    if( current_prefix.size() == counts[dim] )
    {
        shellings_found.push_back( current_prefix );
        LOG << "found" << nl;
        return;
    }

    if( current_prefix.size() == 1 ) 
    {
        current_prefix[0].second = 0;
    }

    // TODO:
    // re-eval the ranks of the remaining_nodes and sort them 
    for( auto& node : remaining_nodes )
    {
        
        const auto faces = mesh.getsubsimplices( dim, dim-1, node.first ).getvalues();
        
        for( int index_f = 0; index_f < faces.size(); index_f++ )
        {
            const int face = faces[index_f];
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( contains( parents, node.first ) );

            if( parents.size() == 1 ) continue;

            for( const auto parent : parents ){
                
                if( parent == node.first ) continue;
                
                for( auto c : current_prefix ) if( c.first == parent ) node.second = minimum( node.second, c.second + 1 );
                
            }
        }
    }

    std::sort( remaining_nodes.begin(), remaining_nodes.end(), 
        []( std::pair<int,Float> x1, std::pair<int,Float> x2 )
        {
            return x1.second < x2.second;
        }
    );


    // There are still volumes to be added
    // Run over all volumes and check whether they can be added 
    for( auto node : remaining_nodes )
    {
        
        assert( not any_of( current_prefix.begin(), current_prefix.end(), [=]( const std::pair<int,Float>& ranked_node ){ return node.first == ranked_node.first; } ) );
        
        
        
        
        {
            for( auto t : current_prefix ) LOG << t.first << " (" << t.second << ")" << space;
            LOG << node.first << " (" << node.second << ")" << nl;    
        }
        

        // list the faces adjacent to the node
        // and check whether it is connected to one of th prefix nodes  
        
        const auto faces_of_node = mesh.getsubsimplices( dim, dim-1, node.first ).getvalues();
        assert( faces_of_node.size() == dim+1 );
        std::vector<bool> face_is_connected( dim+1, false );

        for( int index_f = 0; index_f < faces_of_node.size(); index_f++ )
        {
            const int face = faces_of_node[index_f];
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( 1 <= parents.size() and parents.size() <= 2 );

            if( parents.size() == 1 ) continue;

            assert( parents[0] == node.first or parents[1] == node.first );
            assert( parents[0] != node.first or parents[1] != node.first );

            if( parents[0] != node.first and any_of( current_prefix.begin(), current_prefix.end(), [=]( const std::pair<int,Float>& ranked_node ){ return parents[0] == ranked_node.first; } ) ) face_is_connected[index_f] = true;
            if( parents[1] != node.first and any_of( current_prefix.begin(), current_prefix.end(), [=]( const std::pair<int,Float>& ranked_node ){ return parents[1] == ranked_node.first; } ) ) face_is_connected[index_f] = true;
        }

        
        // OPTIMIZATION
        // if current_prefix is not empty and volume node has no common face with the volumes in current_prefix, then skip 
        
        if( current_prefix.size() != 0 )
        if( not std::any_of( face_is_connected.begin(), face_is_connected.end(), [](bool b) { return b; } ) )
            continue;


        // for each proper subsimplex index_f of the volume s
        // if index_f is a subsimplex of one of the previous volumes 
        // then it must be contained in one of the connected faces 

        bool is_compatible = true;

        // List all proper subsimplices of the node, excluding itself and its faces 
        for( int d = 0; d <= dim-2; d++ )
        {
            if( not is_compatible ) break;
            
            const auto subsimplices_of_node = mesh.getsubsimplices( dim, d, node.first ).getvalues();

            // For each of those proper subsimplices of node ...
            for( const auto sub : subsimplices_of_node )
            {
                
                if( not is_compatible ) break;

                // check whether it has got a supersimplex within the prefix
                const auto supersimplices_of_sub = mesh.getsupersimplices( dim, d, sub );

                bool sub_of_prefix = std::any_of( 
                    supersimplices_of_sub.begin(), supersimplices_of_sub.end(),
                    [&]( int val ) { 
                        
                        return std::any_of( current_prefix.begin(), current_prefix.end(), 
                            [=]( std::pair<int,Float> c ){ 
                                return c.first == val; 
                            }
                        );

                    }
                );

                // if not, there is nothing to check 
                if( not sub_of_prefix ) continue;

                // otherwise, there exist some nodes in the prefix,
                // and we check whether the subsimplex contained within one of the active faces 
                bool within_connected_face = false;

                const auto faces_containing_sub = mesh.getsupersimplices( dim-1, d, sub );

                for( int index_f = 0; index_f < faces_of_node.size(); index_f++ ){
                    if( not face_is_connected[index_f] ) continue;
                    if( contains( faces_containing_sub, faces_of_node[index_f] ) )
                        within_connected_face = true;
                }

                // if sub is not within a connected face, then node is not compatible 
                if( not within_connected_face )
                    is_compatible = false;

            }
        }

        // We have determined whether node is compatible. If not, then continue with the next one

        // is_compatible = true; // TODO

        if( not is_compatible ) {
            LOG << "reject" << nl;
            continue;
        }

        
        
        // What is the dimension of the subsimplex around which the interface is made?

        int k = dim - std::count_if( face_is_connected.begin(), face_is_connected.end(), [](bool b){return b;} );

        if( k < dim-1 )
        {
            // get the subsimplices of the proposed node of dimension k
            const auto k_subsimplices_of_node = mesh.getsubsimplices( dim, k, node.first ).getvalues();

            // find the node which is contained in all connected faces
            int index_of_sub = mesh.nullindex;
            for( int i = 0; i < k_subsimplices_of_node.size(); i++ ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_subsimplices_of_node[i] ) )
                        is_match = false;
                
                if( is_match ) index_of_sub = i;
            }

            assert( index_of_sub != mesh.nullindex );
            for( int i = 0; i < k_subsimplices_of_node.size(); i++ ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_subsimplices_of_node[i] ) )
                        is_match = false;
                if( is_match ) assert( index_of_sub == i );
            }

            // we have found the common subsimplex 
            const int common_subsimplex = k_subsimplices_of_node[index_of_sub];

            // check that all of its parents are here
            const auto& parents = mesh.getsupersimplices( dim, k, common_subsimplex );
            std::vector<bool> parent_is_already_here( parents.size(), false );

            for( int p = 0; p < parents.size(); p++ )
            {
                parent_is_already_here[p] 
                =
                ( parents[p] == node.first ) 
                or
                std::any_of( current_prefix.begin(), current_prefix.end(), 
                    [&]( const std::pair<int,Float>& item ){ return item.first == parents[p]; }
                );
            }

            bool all_parents_are_here = std::all_of( parent_is_already_here.begin(), parent_is_already_here.end(), [](bool b){return b;} );

            assert( all_parents_are_here );

        }

        

        
        
        
        
        
        
        
        // node is compatible with the prefix, then recursion 

        {
            // CHECK that current and remaining nodes are disjoint
            for( const auto& node1 : current_prefix ) 
            for( const auto& node2 : remaining_nodes ) 
                assert( node1.first != node2.first );

            auto next_prefix = current_prefix;
            next_prefix.reserve( counts[dim] );
            next_prefix.push_back( node );

            auto next_remaining_nodes = remaining_nodes;
            auto it = std::find( next_remaining_nodes.begin(), next_remaining_nodes.end(), node );
            // auto it = std::find_if( next_remaining_nodes.begin(), next_remaining_nodes.end(), [=]( std::pair<int,Float> p ) { return p.first == node.first; }  );
            assert( it != next_remaining_nodes.end() );
            next_remaining_nodes.erase( it );

            // CHECK that current and remaining nodes are disjoint
            for( const auto& node1 : next_prefix ) 
            for( const auto& node2 : next_remaining_nodes ) 
                assert( node1.first != node2.first );


            generate_ranked_shelling( mesh, shellings_found, next_prefix, next_remaining_nodes );
        }

    }
    
}









Float estimate_shelling_quality( 
    const Mesh& mesh,
    std::vector<int>
){
    // check input consistency
    const int dim = mesh.getinnerdimension();

}













#endif // SHELLING_LISTER_HPP
