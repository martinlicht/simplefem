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
        // and check whether it is connected to one of the prefix nodes  
        
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
    if( true and shellings_found.size() >= 1 )
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
                assert( current_prefix.size() > 0 );
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

        assert( k < dim or current_prefix.size() > 0 );

        if( k < dim-1 )
        {
            // get the subsimplices of the proposed node of dimension k
            const auto k_subsimplices_of_node = mesh.getsubsimplices( dim, k, node.first ).getvalues();

            // find the node which is contained in all connected faces
            
            int common_subsimplex = mesh.nullindex;
            
            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
                        is_match = false;

                assert( not is_match or common_subsimplex == mesh.nullindex );
                
                if( is_match ) common_subsimplex = k_sub;
            }

            assert( common_subsimplex != mesh.nullindex );

            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
                        is_match = false;

                if( is_match ) assert( common_subsimplex == k_sub );
            }

            // we have found the common subsimplex 
            
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
    std::vector<int> shelling,
    int form_degree
){
    // check input consistency
    
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );
    assert( 0 <= form_degree and form_degree <= dim );

    const auto counts = mesh.count_simplices();

    assert( shelling.size() == counts[dim] );

    // estimate the algebraic condition number of these guys
    
    std::vector<Float> aspect_condition_number(    counts[dim], notanumber );
    std::vector<Float> trafo_singular_max( counts[dim], notanumber );
    std::vector<Float> trafo_singular_min( counts[dim], notanumber );
    std::vector<Float> algebraic_condition_number( counts[dim], notanumber );
    
    for( int i = 0; i < counts[dim]; i++ ) 
    {
        Float diameter = mesh.getDiameter(dim,i);
        
        Float volume = mesh.getMeasure(dim,i);

        const auto faces = mesh.getsubsimplices(dim,dim-1,i).getvalues();

        std::vector<Float> heights( dim+1 );
        for( int f = 0; f < faces.size(); f++ ) 
            heights[f] = dim * volume / mesh.getMeasure( dim-1, faces[f] );

        Float height = *min_element( heights.begin(), heights.end() );

        aspect_condition_number[i] = diameter / height;

        const auto& trafo    = mesh.getTransformationJacobian(dim,i);
        const auto& invtrafo = Inverse(trafo);

        trafo_singular_max[i] =         trafo.operator_norm_estimate();
        trafo_singular_min[i] = 1. / invtrafo.operator_norm_estimate();

        algebraic_condition_number[i] = trafo_singular_max[i] / trafo_singular_min[i];
    }
    
    LOG << "Max aspect condition number:    " << *std::max_element(    aspect_condition_number.begin(),    aspect_condition_number.end() ) << nl;
    LOG << "Max algebraic condition number: " << *std::max_element( algebraic_condition_number.begin(), algebraic_condition_number.end() ) << nl;

    // run over all the n-simplices and compare their diameters
    // Lazy estimate 
    
    Float max_diameter_ratio = 0.;
    for( int e1 = 0; e1 < counts[1]; e1++ ) 
    for( int e2 = 0; e2 < counts[1]; e2++ ) 
    {
        Float diam1 = mesh.getDiameter(1,e1);
        Float diam2 = mesh.getDiameter(1,e2);
        max_diameter_ratio = maximum( max_diameter_ratio, diam1/diam2 );
    }

    LOG << "Max diameter ratio: " << max_diameter_ratio << nl;

    

    

    std::vector<Float> C5(      shelling.size(), notanumber );
    std::vector<Float> C5prime( shelling.size(), notanumber );
    std::vector<Float> C5det(   shelling.size(), notanumber );
    
    std::vector<Float> C6(      shelling.size(), notanumber );
    std::vector<Float> C6prime( shelling.size(), notanumber );
    std::vector<Float> C6det(   shelling.size(), notanumber );
    
    std::vector<Float> C7(      shelling.size(), notanumber );
    std::vector<Float> C7prime( shelling.size(), notanumber );
    std::vector<Float> C7det(   shelling.size(), notanumber );
    
    std::vector<Float> C8(      shelling.size(), notanumber );
    std::vector<Float> C8prime( shelling.size(), notanumber );
    std::vector<Float> C8det(   shelling.size(), notanumber );



    std::vector<FloatVector> coefficient_table( counts[dim], FloatVector( counts[dim], 0. ) );
    
    // start computing the relevant data on each piece 
    for( int i = 0; i < shelling.size(); i++ )
    {
        
        int current_node = shelling[i];
        assert( 0 <= current_node and current_node < counts[dim] );
        
        // list the faces adjacent to the i-th node among the previous nodes 
        
        const auto faces_of_node = mesh.getsubsimplices( dim, dim-1, current_node ).getvalues();
        std::vector<bool> face_is_connected( dim+1, false );
        assert( faces_of_node.size() == face_is_connected.size() );

        for( int index_f = 0; index_f < faces_of_node.size(); index_f++ )
        {
            const int face = faces_of_node[index_f];
            
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( 1 <= parents.size() and parents.size() <= 2 );

            if( parents.size() == 1 ) { 
                assert( mesh.is_subsimplex(dim,dim-1,parents[0],face) ); 
                continue; 
            }

            assert( parents[0] == current_node or parents[1] == current_node );
            assert( parents[0] != current_node or parents[1] != current_node );

            for( int j = 0; j < i; j++ ) 
                face_is_connected[index_f] = face_is_connected[index_f] or ( parents[0] == shelling[j] or parents[1] == shelling[j] );
        }


        // What is the dimension of the subsimplex around which the interface is made?
        
        int k = dim - std::count_if( face_is_connected.begin(), face_is_connected.end(), [](bool b){return b;} );
        Assert( i == 0 or k < dim, i, k );


        // Determine the common subsimplex
        int common_subsimplex = mesh.nullindex; 

        if( i == 0 ) {

            assert( k == dim );

            common_subsimplex = shelling[0];

        } else {
            
            LOGPRINTF( "Subsimplex of dimension %i shared by %i-th simplex, which is %i\n", k, i, current_node );
            
            // get the subsimplices of the proposed node of dimension k
            const auto k_subsimplices_of_node = mesh.getsubsimplices( dim, k, current_node ).getvalues();

            // find the node which is contained in all connected faces
            
            assert( common_subsimplex == mesh.nullindex );
            
            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
                        is_match = false;

                assert( not is_match or common_subsimplex == mesh.nullindex );
                
                if( is_match ) common_subsimplex = k_sub;
            }

            assert( common_subsimplex != mesh.nullindex );

            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
                        is_match = false;

                assert( not is_match or common_subsimplex == k_sub );
            }

            

            // check that all of its parents are here
            const auto& parents = mesh.getsupersimplices( dim, k, common_subsimplex );
            std::vector<bool> parent_is_already_here( parents.size(), false );

            
            LOG << "Parents ";
            for( int p = 0; p < parents.size(); p++ ) LOG << space << parents[p];
            LOG << nl;
            LOG << "Shelling ";
            for( int j = 0; j <= i; j++ ) LOG << shelling[j] << tab;
            LOG << nl;
                
            
            for( int p = 0; p < parents.size(); p++ )
            {
                if ( parents[p] == shelling[i] ) parent_is_already_here[p] = true;

                for( int j = 0; j < i; j++ ) { 
                    if ( parents[p] == shelling[j] ) parent_is_already_here[p] = true;
                }
            }

            bool all_parents_are_here = std::all_of( parent_is_already_here.begin(), parent_is_already_here.end(), [](bool b){return b;} );

            assert( all_parents_are_here );

        }



        // We can now compute all the relevant quantities 
        {
            const   int n = dim;

            const Float B = sqrt( 1. + square( (n-k)*aspect_condition_number[i] ) ) - (n-k)*aspect_condition_number[i];

            const Float Ctheta = max_diameter_ratio;

            const Float kappa = aspect_condition_number[i];

            LOG << i << ":\t" << k << ":\t" << B << space << Ctheta << space << kappa << nl;


            const Float Psi_estimate1 = ( 1. + (k+1) * (1+B) * kappa );
            const Float Psi_estimate = sqrt( 1. + square( (k+1) * (1+B) * kappa ) / 4 ) + ( (k+1) * (1+B) * kappa ) / 2;
            // LOG << Psi_estimate1 << space << Psi_estimate << nl;
            

            
            C5[i]      = (k+1) * Ctheta * kappa * (B*B) * Psi_estimate;
            
            C5prime[i] = (k+1) * B * Ctheta * kappa * B;

            C5det[i]   = power_numerical( (k+1) * B * Ctheta * kappa, n );

            
            C6[i]      = (k+1) * Ctheta * kappa * B * Psi_estimate;
            
            C6prime[i] = (k+1) * B * Ctheta * kappa;

            C6det[i]   = power_numerical( (k+1) * Ctheta * kappa, n );

            
            C7[i]      = ( 1. + 3./2. * (k+1) * kappa ) * (k+1) * kappa * B * Psi_estimate;
            
            C7prime[i] = ( 1. + 3./2. * (k+1) * kappa ) * (k+1) * B * Ctheta * kappa;

            C7det[i]   = 0.5 * power_numerical( (k+1) * B * Ctheta * kappa, n );


            C8[i]      = ( 2. + 3. * (k+1) * kappa ) * (k+1) * kappa * Psi_estimate;
            
            C8prime[i] = ( 2. + 3. * (k+1) * kappa ) * (k+1) * Ctheta * kappa;

            C8det[i]   = 2.0 * power_numerical( (k+1) * Ctheta * kappa, n );
        }

        // Having computed the coefficients, let us now compute the coefficient table

        // fill in values here 
        {
            Float PF = 1. / Constants::pi ; // * 2. * power_numerical( algebraic_condition_number[i], form_degree+1 ) / trafo_singular_min[i];

            Float A = PF;
            Float B = PF * C5[i] * power_numerical( C5prime[i], form_degree   ) * sqrt( C6det[i] );
            Float C =      C5[i] * power_numerical( C5prime[i], form_degree-1 ) * sqrt( C6det[i] );

            // obtain all previous indices of the common subsimplex of dimension n-k
            const auto& relevant_volumes = mesh.getsupersimplices( dim, dim-1, common_subsimplex );

            coefficient_table[i][i] = A;

            for( int j = 0; j < i; j++ )
            {
                bool is_relevant = std::find( relevant_volumes.begin(), relevant_volumes.end(), shelling[j] ) != relevant_volumes.end();

                if( not is_relevant ) continue;

                coefficient_table[i][j] += B;

                coefficient_table[i] += C * coefficient_table[j];
            }

        }


        
        
    }

    // Compute the second estimate 

    Float estimate1 = 0.;

    // skip the first one
    for( int i = 1; i < shelling.size(); i++ )
    for( int j = 1; j < shelling.size(); j++ )
    {
        estimate1 += square( coefficient_table[i][j] );
    }

    LOG << "First estimate: " << sqrt( estimate1 ) << nl;


    Float estimate2 = 1.;

    // skip the first one
    for( int i = 1; i < shelling.size(); i++ )
    {
        estimate2 *= C8[i] * power_numerical( C8prime[i], form_degree ) * C7[i] * power_numerical( C7prime[i], form_degree-1 ) * sqrt( C8det[i] * C7det[i] );
        LOG << i << tab << estimate2 << nl;
    }

    LOG << "Second estimate: " << estimate2 << nl;

    return estimate1;

}













#endif // SHELLING_LISTER_HPP
