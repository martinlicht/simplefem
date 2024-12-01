#ifndef SHELLING_LISTER_HPP
#define SHELLING_LISTER_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <utility> // For std::pair

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "mesh.hpp"




void generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list,
    std::vector<std::vector<int>>& shellings_found,
    std::vector<int> current_prefix,
    std::vector<int> remaining_nodes
){
    
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );


    // count all the simplices of the mesh 
    const auto counts = mesh.count_simplices();


    // TODO:
    // if enough shellings have been found already, then return 
    if( shellings_found.size() >= 1 )
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
        
        const auto faces = mesh.getsubsimplices( dim, dim-1, node ).getvalues();
        assert( faces.size() == dim+1 );
        std::vector<bool> face_is_connected( dim+1, false );

        for( int index_f = 0; index_f < faces.size(); index_f++ )
        {
            const int face = faces[index_f];
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( 1 <= parents.size() and parents.size() <= 2 );

            if( parents.size() == 1 ) continue;

            assert( parents[0] == node or parents[1] == node );

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
            const auto face = faces[index_f];
            if( current_prefix.size() == 0 )
            if( acceptable_face_list[face] and face_is_connected[index_f] ) 
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

                for( int index_f = 0; index_f < faces.size(); index_f++ ){
                    if( not face_is_connected[index_f] ) continue;
                    if( contains( faces_containing_sub, faces[index_f] ) )
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
            next_remaining_nodes.erase( it );

            generateShellings( mesh, acceptable_face_list, shellings_found, next_prefix, next_remaining_nodes );
        }

    }
    

    
}






std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh,
    const std::vector<bool>& acceptable_face_list
){
    const int dim = mesh.getinnerdimension();

    std::vector<std::vector<int>> shellings_found;
    
    std::vector<int> current_prefix;

    std::vector<int> remaining_nodes( mesh.count_simplices( dim ) );
    for( int node = 0; node < remaining_nodes.size(); node++ ) remaining_nodes[node] = node;

    generateShellings( mesh, acceptable_face_list, shellings_found, current_prefix, remaining_nodes );
    
    return shellings_found;
}



std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh
){
    const int dim = mesh.getinnerdimension();

    std::vector<bool> acceptable_face_list( mesh.count_simplices(dim) );
    
    for( int face = 0; face < mesh.count_simplices(dim-1); face++ ) 
    {
        const auto& parents = mesh.getsupersimplices(dim,dim-1,face);
        acceptable_face_list[face] = ( parents.size() != 1 );
    }

    return generateShellings( mesh, acceptable_face_list );
}











#endif // SHELLING_LISTER_HPP
