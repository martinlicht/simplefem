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
    const std::vector<std::vector<int>>& volume_acceptable_face_list,
    std::vector<std::vector<int>>& shellings_found,
    std::vector<int> current_prefix,
    std::vector<int> remaining_nodes
){
    
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const auto counts = mesh.count_simplices();

    // if the current prefix is already containing all volumes, then we are done 
    if( current_prefix.size() == counts[dim] )
    {
        shellings_found.push_back( current_prefix );
        LOG << "found" << nl;
        return;
    }

    // there are still volumes to be added
    // run over all volumes and check whether they can be added 
    // for( int s = 0; s < counts[dim]; s++ )
    for( int s : remaining_nodes )
    {
        
        // if already contained within prefix, then skip 
        if( contains( current_prefix, s ) ) 
            continue;
        
        
        {
            for( int t : current_prefix ) LOG << t << space;
            LOG << s << nl;    
        }
        

        // list the faces and check whether it is connected to one of the previous simplices 
        
        const auto faces = mesh.getsubsimplices(dim,dim-1,s).getvalues();
        assert( faces.size() == dim+1 );
        std::vector<bool> face_is_connected( dim+1, false );

        for( int f = 0; f < faces.size(); f++ )
        {
            const int face = faces[f];
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( contains( parents, s ) );

            for( const auto parent : parents ){
                if( parent == s ) continue;
                if( contains( current_prefix, parent ) ) face_is_connected[f] = true;
            }
        }

        
        // OPTIMIZATION
        // if current_prefix is not empty and volume s has no common face with the volumes in current_prefix, then skip 
        
        if( current_prefix.size() != 0 )
        if( not std::any_of( face_is_connected.begin(), face_is_connected.end(), [](bool val) { return val; } ) )
            continue;


        // if the connection to the previous volumes is not via an acceptable face, then skip
        bool acceptable_connection = false;
        for( int f = 0; f <= dim; f++ ){
            if( current_prefix.size() > 0 and not face_is_connected[f] ) continue;
            if( contains( volume_acceptable_face_list[s], faces[f] ) )
                acceptable_connection = true;
        }

        if( not acceptable_connection ) {
            LOG << "unacceptable" << nl;
            continue;
        }
        
        
        // for each proper subsimplex f of the volume s
        // if f is a subsimplex of one of the previous volumes 
        // then it must be contained in one of the connected faces 

        bool is_compatible = true;

        // List all proper subsimplices of s, excluding s and its faces 
        for( int d = 0; d <= dim-2; d++ )
        {
            if( not is_compatible ) break;
            
            const auto subsimplices_of_s = mesh.getsubsimplices( dim, d, s ).getvalues();

            // For each of those proper subsimplices of s ...
            for( const auto sub : subsimplices_of_s )
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

                for( int f = 0; f < faces.size(); f++ ){
                    if( not face_is_connected[f] ) continue;
                    if( contains( faces_containing_sub, faces[f] ) )
                        within_connected_face = true;
                }

                // if sub is not within a connected face, then s is not compatible 
                if( not within_connected_face )
                    is_compatible = false;

            }
        }

        // We have determined whether s is compatible. If not, then continue with the next one

        // is_compatible = true; // TODO

        if( not is_compatible ) {
            LOG << "reject" << nl;
            continue;
        }

        // s is compatible with the prefix, then recursion 

        {
            auto next_prefix = current_prefix;
            next_prefix.reserve( counts[dim] );
            next_prefix.push_back( s );

            auto next_remaining_nodes = remaining_nodes;
            auto it = std::find( next_remaining_nodes.begin(), next_remaining_nodes.end(), s );
            next_remaining_nodes.erase( it );

            generateShellings( mesh, volume_acceptable_face_list, shellings_found, next_prefix, next_remaining_nodes );
        }

    }
    

    
}






std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh,
    const std::vector<std::vector<int>>& volume_acceptable_face_list
){
    const int dim = mesh.getinnerdimension();

    std::vector<std::vector<int>> shellings_found;
    
    std::vector<int> current_prefix;

    std::vector<int> remaining_nodes( mesh.count_simplices( dim ) );
    for( int s = 0; s < remaining_nodes.size(); s++ ) remaining_nodes[s] = s;

    generateShellings( mesh, volume_acceptable_face_list, shellings_found, current_prefix, remaining_nodes );
    
    return shellings_found;
}



std::vector<std::vector<int>> generateShellings(
    const Mesh& mesh
){
    const int dim = mesh.getinnerdimension();

    // create the list of faces of each volume that is acceptable 
    std::vector<std::vector<int>> volume_acceptable_face_list( mesh.count_simplices(dim) );
    
    for( int s = 0; s < mesh.count_simplices(dim); s++ ) 
    {
        const auto faces = mesh.getsubsimplices(dim,dim-1,s).getvalues();
        
        for( auto face : faces ) 
        {
            const auto& parents = mesh.getsupersimplices(dim,dim-1,face);
            if( parents.size() == 1 ) continue;

            assert( 0 < parents.size() && parents.size() <= 2 );
            
            volume_acceptable_face_list[s].push_back(face);
        }
    }

    return generateShellings( mesh, volume_acceptable_face_list );
}











#endif // SHELLING_LISTER_HPP
