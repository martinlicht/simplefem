#ifndef SHELLING_LISTER_2_HPP
#define SHELLING_LISTER_2_HPP

#include <algorithm>
#include <vector>
#include <queue>
#include <utility> // For std::pair

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "mesh.hpp"





struct mesh_information_for_shelling
{
    Float max_diameter_ratio = notanumber;
    std::vector<Float> diameters;
    std::vector<Float> volumes;
    std::vector<Float> heights;
    
    std::vector<Float> trafo_singular_max;
    std::vector<Float> trafo_singular_min;
    
    std::vector<Float> aspect_condition_number;
    std::vector<Float> algebraic_condition_number;
    
    std::vector<std::vector<std::vector<Float>>> C5;
    std::vector<std::vector<std::vector<Float>>> C6;
    std::vector<std::vector<std::vector<Float>>> C7;
    std::vector<std::vector<std::vector<Float>>> C8;
    
    // Constructor
    mesh_information_for_shelling( const Mesh& mesh )
    {
        
        // check input consistency
    
        const int dim = mesh.getinnerdimension();
        assert( dim >= 1 );

        const auto counts = mesh.count_simplices();

        // Initialize vectors with required sizes and fill them with notanumber
        
        diameters.resize(counts[dim], notanumber);
        volumes.resize(counts[dim], notanumber);
        heights.resize(counts[dim], notanumber);
        
        trafo_singular_max.resize(counts[dim], notanumber);
        trafo_singular_min.resize(counts[dim], notanumber);
        
        aspect_condition_number.resize(counts[dim], notanumber);
        algebraic_condition_number.resize(counts[dim], notanumber);
        
        C5.resize( counts[dim]+1, std::vector<std::vector<Float>>( dim+1, std::vector<Float>( dim+1, notanumber ) ) );
        C6.resize( counts[dim]+1, std::vector<std::vector<Float>>( dim+1, std::vector<Float>( dim+1, notanumber ) ) );
        C7.resize( counts[dim]+1, std::vector<std::vector<Float>>( dim+1, std::vector<Float>( dim+1, notanumber ) ) );
        C8.resize( counts[dim]+1, std::vector<std::vector<Float>>( dim+1, std::vector<Float>( dim+1, notanumber ) ) );
        
        // estimate several geometric properties of each n-simplex
        // - diameter
        // - volume
        // - aspect condition number
        // - algebraic condition number
    
        for( int i = 0; i < counts[dim]; i++ ) 
        {
            diameters[i] = mesh.getDiameter(dim,i);
            
            volumes[i] = mesh.getMeasure(dim,i);

            const auto faces = mesh.getsubsimplices(dim,dim-1,i).getvalues();

            std::vector<Float> heights( dim+1 );
            for( int f = 0; f < faces.size(); f++ ) 
                heights[f] = dim * volumes[i] / mesh.getMeasure( dim-1, faces[f] );

            Float min_height = *min_element( heights.begin(), heights.end() );

            aspect_condition_number[i] = diameters[i] / min_height;

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


        // Compute the coefficients 
        for( int i = 0; i < counts[dim]; i++ ) 
        for( int l = 0; l <= dim; l++ )
        {
            const   int n = dim;

            const Float Ctheta = max_diameter_ratio;

            const Float kappa = aspect_condition_number[i];

            const Float mu_l_part = (n-l) * kappa;
            
            const Float mu_l = sqrt( 1. + square( mu_l_part ) ) + mu_l_part;

            const Float nu_l_part = (l+1) * (1+mu_l) * kappa / 2.;
            
            const Float nu_l = sqrt( 1. + square( nu_l_part ) ) + nu_l_part;

            for( int k = 0; k <= dim; k++ )
            {
                
                C5[i][l][k]   = mu_l * nu_l * power_numerical( mu_l * kappa * Ctheta * (l+1), k - n/2. );
            
                C6[i][l][k]   = mu_l * nu_l * power_numerical(        kappa * Ctheta * (l+1), k - n/2. );

                const Float temp = 1. + 3./2. * (l+1) * kappa;
            
                C7[i][l][k]   = temp * nu_l * power_numerical( mu_l * kappa * Ctheta * (l+1), k - n/2. ) * power_numerical(2,n/2.);
            
                C8[i][l][k]   = temp * nu_l * power_numerical(        kappa * Ctheta * (l+1), k - n/2. ) * power_numerical(2,n/2.);

        }

        for( auto a : C5 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) );
        for( auto a : C6 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) );
        for( auto a : C7 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) );
        for( auto a : C8 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) );
        
        for( auto a : aspect_condition_number ) assert( not std::isnan(a) );

        for( auto a : algebraic_condition_number ) assert( not std::isnan(a) );
        
    }
};





std::vector<std::vector<int>> generate_shellings( 
    const Mesh& mesh,
    int form_degree
);

std::vector<std::vector<int>> generate_shellings(
    const Mesh& mesh,
    int form_degree,
    const mesh_information_for_shelling& info,
    std::vector<std::vector<int>>& shellings_found,
    const std::vector<int>& current_prefix,
    const std::vector<FloatVector>& coefficient_table
);

std::vector<std::vector<int>> generate_shellings( 
    const Mesh& mesh,
    int form_degree
){
    
    const auto info = mesh_information_for_shelling( mesh );

    std::vector<std::vector<int>> shellings_found;

    std::vector<int> current_prefix;

    std::vector<FloatVector> coefficient_table;

    auto result = generate_shellings( mesh, form_degree, info, shellings_found, current_prefix, coefficient_table );

    return result;
}


std::vector<std::vector<int>> generate_shellings(
    const Mesh& mesh,
    int form_degree,
    const mesh_information_for_shelling& info,
    std::vector<std::vector<int>>& shellings_found,
    const std::vector<int>& current_prefix,
    const std::vector<FloatVector>& coefficient_table
){

    const int dim = mesh.getinnerdimension();

    const auto counts = mesh.count_simplices();

    
    // check that input mesh is reasonable 
    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );


    // count all the simplices of the mesh 
    const auto counts = mesh.count_simplices();


    // if enough shellings have been found already, then return 
    if( false and shellings_found.size() >= 1 )
        return;


    // if the current prefix is already containing all volumes, 
    // then a shelling has been found and can be added to the list
    // we can return from this place, there is nothing to do here 
    if( current_prefix.size() == counts[dim] )
    {
        shellings_found.push_back( current_prefix );
        LOG << "found" << nl;
        return;
    }
    assert( current_prefix.size() < counts[dim]               );
    assert( current_prefix.size() == coefficient_table.size() );



    // collect all the nodes not yet in the shelling
    
    std::vector<int> remaining_nodes;
    for( int s = 0; s < counts[dim]; s++ )
    {
        bool already_in_prefix = std::any_of( current_prefix.begin(), current_prefix.end(), [=]( const int node ){ return node == s; } );

        if( current_prefix.size() == 0 ) assert( not already_in_prefix );

        if( not already_in_prefix ) remaining_nodes.push_back( s );
    }

    assert( remaining_nodes.size() > 0                                    );
    assert( remaining_nodes.size() + current_prefix.size() == counts[dim] );
    for( int s : remaining_nodes ) for( int p : current_prefix ) assert( s != p );
    
    
    // for each unused node, collect the active faces 

    std::vector<int>               how_many_connected_faces( remaining_nodes.size(), -1                         );
    std::vector<std::vector<bool>> face_is_connected       ( remaining_nodes.size(), std::vector<bool>( dim+1 ) );

    for( int i = 0; i < remaining_nodes.size(); i++ )
    {
        int node = remaining_nodes[i];
        
        assert( 0 <= node and node < counts[dim] );
        
        const auto faces_of_node = mesh.getsubsimplices( dim, dim-1, node ).getvalues();

        assert( faces_of_node.size() == face_is_connected[i].size() );

        how_many_connected_faces[i] = 0;

        for( int index_f = 0; index_f < faces_of_node.size(); index_f++ )
        {
            const int face = faces_of_node[index_f];
            
            const auto parents = mesh.getsupersimplices( dim, dim-1, face );
            
            assert( 1 <= parents.size() and parents.size() <= 2 );

            if( parents.size() == 1 ) { 
                assert( mesh.is_subsimplex(dim,dim-1,parents[0],face) ); 
                face_is_connected[i][index_f] = false;
                continue; 
            }

            assert( parents[0] == node or parents[1] == node );
            assert( parents[0] != node or parents[1] != node );

            for( int p : current_prefix ) 
                face_is_connected[i][index_f] = face_is_connected[i][index_f] or ( parents[0] == p or parents[1] == p );

            if( face_is_connected[i][index_f] ) how_many_connected_faces[i]++;
        }

        assert( how_many_connected_faces[i] == std::count_if( face_is_connected[i].begin(), face_is_connected[i].end(), [=](bool b){return b;} ) );
    
    }


    // unless we are at the start, delete all the non-reachable nodes, that is, those without active faces 
    if( current_prefix.size() > 0 )
    {
        int i = 0;
        while( i < remaining_nodes.size() )
        {
            if( how_many_connected_faces[i] > 0 ) {
                i++;
            } else {
                remaining_nodes.erase(          remaining_nodes.begin() + i          );
                face_is_connected.erase(        face_is_connected.begin() + i        );
                how_many_connected_faces.erase( how_many_connected_faces.begin() + i );
            }
        }

        assert( remaining_nodes.size() > 0 );
        assert( remaining_nodes.size() == face_is_connected.size()        );
        assert( remaining_nodes.size() == how_many_connected_faces.size() );
    }




    
    
    
    // compute the weight associated with the extension along each possible node 

    std::vector<int>   shared_subsimplex_dim( remaining_nodes.size(), -1    );
    std::vector<bool>  shelling_compatible  ( remaining_nodes.size(), true );

    std::vector<Float> weight_for_node_1( remaining_nodes.size(), std::numeric_limits<Float>::infinity() );
    std::vector<Float> weight_for_node_2( remaining_nodes.size(), std::numeric_limits<Float>::infinity() );

    for( int i = 0; i < remaining_nodes.size(); i++ ) 
        assert( std::isfinite( weight_for_node_1[i] ) == std::isfinite( weight_for_node_1[i] ) );



    // start computing the relevant data on each node  
    for( int i = 0; i < remaining_nodes.size(); i++ ) 
    {
        
        int current_node = remaining_nodes[i];
        assert( 0 <= current_node and current_node < counts[dim] );
        
        // What is the dimension of the subsimplex around which the interface is made?
        
        int k = dim - how_many_connected_faces[i];
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
            
            const auto faces_of_node = mesh.getsubsimplices( dim, dim-1, current_node ).getvalues();

            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[i][f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
                        is_match = false;

                assert( not is_match or common_subsimplex == mesh.nullindex );
                
                if( is_match ) common_subsimplex = k_sub;
            }

            assert( common_subsimplex != mesh.nullindex );

            for( int k_sub : k_subsimplices_of_node ) 
            {
                bool is_match = true;
                for( int f = 0; f < faces_of_node.size() && is_match; f++ )
                    if( face_is_connected[i][f] and not mesh.is_subsimplex( dim-1, k, faces_of_node[f], k_sub ) )
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

    
    
    // delete all the nodes that don't give a shelling 
    if( current_prefix.size() > 0 )
    {
        int i = 0;
        while( i < remaining_nodes.size() )
        {
            if( shelling_compatible[i] ) {
                i++;
            } else {
                remaining_nodes.erase(          remaining_nodes.begin() + i          );
                face_is_connected.erase(        face_is_connected.begin() + i        );
                how_many_connected_faces.erase( how_many_connected_faces.begin() + i );
                shared_subsimplex_dim.erase(    shared_subsimplex_dim.begin() + i    );
                shelling_compatible.erase(      shelling_compatible.begin() + i      );
                weight_for_node_1.erase(        weight_for_node_1.begin() + i        );
                weight_for_node_2.erase(        weight_for_node_2.begin() + i        );
            }
        }

        assert( remaining_nodes.size() > 0 );
        assert( remaining_nodes.size() == face_is_connected.size()        );
        assert( remaining_nodes.size() == how_many_connected_faces.size() );
        assert( remaining_nodes.size() == shared_subsimplex_dim.size()    );
        assert( remaining_nodes.size() == shelling_compatible.size()      );
        assert( remaining_nodes.size() == weight_for_node_1.size()        );
        assert( remaining_nodes.size() == weight_for_node_2.size()        );
    }

    // sort nodes by priority
    // simply some lazy bubble sort
    for( int i = 0; i < remaining_nodes.size(); i++ )
    for( int j = 0; j < remaining_nodes.size(); j++ )
    {
        bool correct_order = weight_for_node_1[i] < weight_for_node_1[j];

        if( correct_order ) continue;

        std::swap( remaining_nodes[i],          remaining_nodes[j]          );
        std::swap( how_many_connected_faces[i], how_many_connected_faces[j] );
        std::swap( face_is_connected[i],        face_is_connected[j]        );
        std::swap( shared_subsimplex_dim[i],    shared_subsimplex_dim[j]    );
        std::swap( shelling_compatible[i],      shelling_compatible[j]      );
        std::swap( weight_for_node_1[i],        weight_for_node_1[j]        );
        std::swap( weight_for_node_2[i],        weight_for_node_2[j]        );
    }
    

    // only reachable nodes are left and they are ordered by priority 

    for( int i = 0; i < remaining_nodes.size(); i++ )
    {
        int node = remaining_nodes[i];
        
        auto next_prefix = current_prefix;
        next_prefix.push_back( node );

        generate_shellings( mesh, form_degree, info, shellings_found, next_prefix, coefficient_table );
    }



    // Compute the second estimate 

    Float estimate1 = 0.;

    // skip the first one
    for( int i = 1; i < counts[dim]; i++ )
    for( int j = 1; j < counts[dim]; j++ )
    {
        estimate1 += square( coefficient_table[i][j] );
    }

    LOG << "First estimate: " << sqrt( estimate1 ) << nl;


    Float estimate2 = 1.;

    // skip the first one
    for( int i = 1; i < counts[dim]; i++ )
    {
        estimate2 *= C8[i] * power_numerical( C8prime[i], form_degree ) * C7[i] * power_numerical( C7prime[i], form_degree-1 ) * sqrt( C8det[i] * C7det[i] );
        LOG << i << tab << estimate2 << nl;
    }

    LOG << "Second estimate: " << estimate2 << nl;

    return estimate1;

}













#endif // SHELLING_LISTER_2_HPP
