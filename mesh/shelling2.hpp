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
    Float min_height_of_vertex = notanumber;
    Float max_volume_ratio   = notanumber;

    std::vector<Float> diameters;
    std::vector<Float> volumes;
    std::vector<Float> heights;
    
    std::vector<Float> trafo_singular_max;
    std::vector<Float> trafo_singular_min;
    
    std::vector<Float> aspect_condition_number;
    std::vector<Float> algebraic_condition_number;

    std::vector<std::vector<FloatVector>> midpoints;            // [dimension][simplex index]
    std::vector<std::vector<Float>>       heights_of_vertex;    // [volume index][local vertex index]
    std::vector<std::vector<FloatVector>> heightvectors_of_vertex;    // [volume index][local vertex index]
    
    std::vector<std::vector<std::vector<Float>>> C5;
    std::vector<std::vector<std::vector<Float>>> C6;
    std::vector<std::vector<std::vector<Float>>> C7;
    std::vector<std::vector<std::vector<Float>>> C8;
    
    // Constructor
    mesh_information_for_shelling( const Mesh& mesh );

};

mesh_information_for_shelling::mesh_information_for_shelling( const Mesh& mesh )
{
        
    // check input consistency

    const int dim = mesh.getinnerdimension();
    assert( dim >= 1 );

    const auto counts = mesh.count_simplices();

    // collect the midpoints of all simplices 
    for( int d = 0; d <= dim; d++ )
    {
        std::vector<FloatVector> midpoints_d;
        midpoints_d.reserve( counts[dim] );

        for( int s = 0; s < counts[d]; s++ ) {
            auto midpoint = mesh.get_midpoint( d, s );
            assert( midpoint.isfinite() );
            midpoints_d.push_back( midpoint );
        }

        midpoints.push_back( midpoints_d );
    }
    assert( midpoints.size() == dim+1 );
    for( int d = 0; d < midpoints.size(); d++ ) assert( midpoints[d].size() == counts[d] );
    
    // Initialize vectors with required sizes and fill them with notanumber
    
    diameters.resize(counts[dim], notanumber);
    volumes.resize(counts[dim], notanumber);
    heights.resize(counts[dim], notanumber);
    
    trafo_singular_max.resize(counts[dim], notanumber);
    trafo_singular_min.resize(counts[dim], notanumber);
    
    aspect_condition_number.resize(counts[dim], notanumber);
    algebraic_condition_number.resize(counts[dim], notanumber);
    
    // estimate several geometric properties of each n-simplex
    // - diameter
    // - volume
    // - aspect condition number
    // - algebraic condition number
    // - minimum height 
    // - heights of all vertices in each simplex 

    min_height_of_vertex = std::numeric_limits<Float>::infinity();

    heights_of_vertex.resize( counts[dim], std::vector<Float>( dim+1, notanumber ) );
    
    for( int i = 0; i < counts[dim]; i++ ) 
    {
        // collect diameter and volume 

        diameters[i] = mesh.getDiameter(dim,i);
        
        volumes[i] = mesh.getMeasure(dim,i);

        // collect the heights 

        const auto faces = mesh.getsubsimplices(dim,dim-1,i).getvalues();

        std::vector<Float> heights( dim+1 );
        for( int f = 0; f < faces.size(); f++ ) 
            heights[f] = dim * volumes[i] / mesh.getMeasure( dim-1, faces[f] );

        Float min_height = *min_element( heights.begin(), heights.end() );

        min_height_of_vertex = minimum( min_height, min_height_of_vertex );

        for( int f = 0; f <= dim; f++ ) {
            int vi = mesh.get_opposite_subsimplex_index( dim, dim-1, i, f );
            heights_of_vertex[i][vi] = heights[f];
        }

        for( int vi = 0; vi <= dim; vi++ ) assert( std::isfinite( heights_of_vertex[i][vi] ) );

        
        // estimate the aspect / algebraic condition number 

        aspect_condition_number[i] = diameters[i] / min_height;

        const auto& trafo    = mesh.getTransformationJacobian(dim,i);
        const auto& invtrafo = Inverse(trafo);

        trafo_singular_max[i] =         trafo.operator_norm_estimate();
        trafo_singular_min[i] = 1. / invtrafo.operator_norm_estimate();

        assert( trafo.isfinite() );
        assert( invtrafo.isfinite() );
        Assert( trafo_singular_max[i] > 0., trafo_singular_max[i] );
        assert( trafo_singular_min[i] > 0. );

        algebraic_condition_number[i] = trafo_singular_max[i] / trafo_singular_min[i];
    }
        
    LOG << "Max aspect condition number:    " << *std::max_element(    aspect_condition_number.begin(),    aspect_condition_number.end() ) << nl;
    LOG << "Max algebraic condition number: " << *std::max_element( algebraic_condition_number.begin(), algebraic_condition_number.end() ) << nl;
    
    
    heightvectors_of_vertex.resize( counts[dim], std::vector<FloatVector>( dim+1, FloatVector( dim, notanumber ) ) );

    for( int i = 0; i < counts[dim]; i++ ) 
    for( int d = 0; d <= dim; d++ )
    {
        heightvectors_of_vertex[i][d] = mesh.getHeightVector( dim, i, d );
    }

    for( auto a : heightvectors_of_vertex ) for( auto b : a ) assert( b.isfinite() );
    
    
    
    // run over all the n-simplices and compare their diameters
    // Lazy estimate 
    
    max_diameter_ratio = 0.;
    
    for( int e1 = 0; e1 < counts[1]; e1++ ) 
    for( int e2 = 0; e2 < counts[1]; e2++ ) 
    {
        Float diam1 = mesh.getDiameter(1,e1);
        Float diam2 = mesh.getDiameter(1,e2);
        max_diameter_ratio = maximum( max_diameter_ratio, diam1/diam2 );
    }

    LOG << "Max diameter ratio: " << max_diameter_ratio << nl;

    Float max_volume_ratio = 0.;
    for( int f = 0; f < counts[dim-1]; f++ )
    {
        const auto parents = mesh.getsupersimplices(dim,dim-1,f);
        if( parents.size() == 1 ) continue;
        assert( parents.size() == 2 );
        const auto p0 = parents[0];
        const auto p1 = parents[1];
        Float v0 = mesh.getMeasure( dim-1, p0 );
        Float v1 = mesh.getMeasure( dim-1, p1 );
        Float v1v0 = v1 / v0;
        Float v0v1 = v0 / v1;
        max_volume_ratio = maximum( max_volume_ratio, v1v0, v0v1 );
    }

    C5.resize( counts[dim], std::vector<std::vector<Float>>( dim, std::vector<Float>( dim+1, notanumber ) ) );
    C6.resize( counts[dim], std::vector<std::vector<Float>>( dim, std::vector<Float>( dim+1, notanumber ) ) );
    C7.resize( counts[dim], std::vector<std::vector<Float>>( dim, std::vector<Float>( dim+1, notanumber ) ) );
    C8.resize( counts[dim], std::vector<std::vector<Float>>( dim, std::vector<Float>( dim+1, notanumber ) ) );

    // Compute the coefficients 
    for( int i = 0; i < counts[dim]; i++ ) 
    for( int l = 0; l < dim; l++ )
    {
        const   int n = dim;

        const Float Ctheta = max_diameter_ratio;

        const Float kappa  = aspect_condition_number[i];

        const Float rho    = minimum( 1., min_height_of_vertex / (l+1.) );

        

        const Float mu_l = 0.5 * sqrt( square(1+rho) + square( square(1+rho) * (n-l) * kappa ) ) + 0.5 * sqrt( square(1-rho) + square( square(1+rho) * (n-l) * kappa ) );
        
        const Float xi_1 = sqrt( square(2*rho+1) + square( (n-l) * kappa ) ) + sqrt( 1. + square( (n-l) * kappa ) );
        
        
        assert( not std::isnan( kappa  ) && kappa  > 0 );
        assert( not std::isnan( Ctheta ) && Ctheta > 0 );
        assert( not std::isnan( mu_l   ) && mu_l   > 0 );
        
        for( int k = 0; k <= dim; k++ )
        {
            Float smax    = mu_l;
            Float sinvmin = mu_l / rho;
            Float sdet    = rho;

            // if(false)
            if( l == dim-1 ) {
                smax = 0.5 * ( sqrt( square( Ctheta * kappa + 1. ) + kappa ) + sqrt( square( Ctheta * kappa - 1. ) + kappa ) );
                sinvmin = smax / max_volume_ratio;
                sdet = max_volume_ratio;
            }
            
            if ( k == 0 ) {
                C5[i][l][k] = power_numerical( sdet, -n/2. );
                C6[i][l][k] = power_numerical( sdet, +n/2. );
            } else if( 0 < k && k < dim ) {
                C5[i][l][k] = smax    * power_numerical( sdet, -n/2. );
                C6[i][l][k] = sinvmin * power_numerical( sdet, +n/2. );
            } else {
                assert( k == dim );
                C5[i][l][k] = power_numerical( sdet, +1.-n/2. );
                C6[i][l][k] = power_numerical( sdet, -1.+n/2. );
            }

            Float tmax    = xi_1 / ( 2 * (1.+rho) );
            Float tinvmin = xi_1 / ( 2 * rho      );
            Float tdet    =  rho / ( 1 + rho );
            
            if ( k == 0 ) {
                C7[i][l][k] = power_numerical( tdet, -n/2. );
                C8[i][l][k] = power_numerical( tdet, +n/2. );
            } else if( 0 < k && k < dim ) {
                C7[i][l][k] = tmax    * power_numerical( tdet, -n/2. );
                C8[i][l][k] = tinvmin * power_numerical( tdet, +n/2. );
            } else {
                assert( k == dim );
                C7[i][l][k] = power_numerical( tdet, +1.-n/2. );
                C8[i][l][k] = power_numerical( tdet, -1.+n/2. );
            }

            Assert( not std::isnan( C5[i][l][k] ), k ); Assert( C5[i][l][k] > 0. ); 
            Assert( not std::isnan( C6[i][l][k] ), k ); Assert( C6[i][l][k] > 0. ); 
            Assert( not std::isnan( C7[i][l][k] ), k ); Assert( C7[i][l][k] > 0. ); 
            Assert( not std::isnan( C8[i][l][k] ), k ); Assert( C8[i][l][k] > 0. ); 
            
        }

    }

    for( auto a : C5 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) && c > 0. );
    for( auto a : C6 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) && c > 0. );
    for( auto a : C7 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) && c > 0. );
    for( auto a : C8 ) for( auto b : a ) for( auto c : b ) assert( not std::isnan(c) && c > 0. );
    
    for( auto a : aspect_condition_number ) assert( not std::isnan(a) );

    for( auto a : algebraic_condition_number ) assert( not std::isnan(a) );

    LOG << "Finished collecting information about the mesh." << nl;
    LOGPRINTF("%e %e %e\n", max_volume_ratio, min_height_of_vertex, max_diameter_ratio );
    
}









struct shelling
: public std::vector<int>
{
    Float weight_reflection;
    Float weight_deformation;

    shelling( const std::vector<int>& indices, Float weight_reflection, Float weight_deformation )
    : std::vector<int>( indices ), weight_reflection( weight_reflection ), weight_deformation( weight_deformation )
    {}
};













std::vector<shelling> generate_shellings2( 
    const Mesh& mesh,
    int form_degree
);

void generate_shellings2(
    const Mesh& mesh,
    int form_degree,
    const mesh_information_for_shelling& info,
    std::vector<shelling>& shellings_found,
    const std::vector<int>& current_prefix,
    const std::vector<FloatVector>& coefficient_table,
    Float multweight
);

std::vector<shelling> generate_shellings2( 
    const Mesh& mesh,
    int form_degree
){
    
    const auto info = mesh_information_for_shelling( mesh );

    std::vector<shelling> shellings_found;

    std::vector<int> current_prefix;

    std::vector<FloatVector> coefficient_table;

    Float multweight = 1.;

    current_prefix.reserve( mesh.count_simplices(mesh.getinnerdimension()) );

    generate_shellings2( mesh, form_degree, info, shellings_found, current_prefix, coefficient_table, multweight );

    return shellings_found;
}


void generate_shellings2(
    const Mesh& mesh,
    int form_degree,
    const mesh_information_for_shelling& info,
    std::vector<shelling>& shellings_found,
    const std::vector<int>& current_prefix,
    const std::vector<FloatVector>& coefficient_table,
    Float multweight
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
    // then a shelling has been found and can be added to the list
    // we can return from this place, there is nothing to do here 
    if( current_prefix.size() == counts[dim] )
    {
        // Compute the estimates  

        Float estimate1 = 0.;

        LOG << "FOUND: ";
        for( int i = 0; i < counts[dim]; i++ ) {
            estimate1 += coefficient_table[i].norm_sq();
            LOGPRINTF( "%i %.4f\t", current_prefix[i], estimate1 );
        }
        estimate1 = sqrt(estimate1);
        LOGPRINTF( "w1=%e w2=%e\n", estimate1, multweight );
        
        shellings_found.push_back( shelling( current_prefix, estimate1, multweight ) );

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
    
    
    // for each as-of-now unused node, collect the active faces 

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

    // unless we are at the start: 
    // delete all the non-reachable nodes, that is, those without active faces 
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

    // LOG << "remaining nodes: "; for( auto node : remaining_nodes ) LOG << node << space; LOG << nl;

    


    
    
    
    // compute the weight associated with the extension along each possible node 

    assert( remaining_nodes.size() > 0 );

    std::vector<int>   shared_subsimplex_dim( remaining_nodes.size(), -1    );
    std::vector<bool>  shelling_compatible  ( remaining_nodes.size(), true );

    std::vector<Float> weight_for_node_reflect( remaining_nodes.size(), std::numeric_limits<Float>::infinity() );
    std::vector<Float> weight_for_node_morphin( remaining_nodes.size(), std::numeric_limits<Float>::infinity() );
    
    std::vector<FloatVector> coefficient_vectors( remaining_nodes.size(), FloatVector( counts[dim], notanumber ) );

    for( int i = 0; i < remaining_nodes.size(); i++ ) 
        assert( std::isfinite( weight_for_node_reflect[i] ) == std::isfinite( weight_for_node_reflect[i] ) );



    // start computing the relevant data on each node  
    for( int i = 0; i < remaining_nodes.size(); i++ ) 
    {
        
        int current_node = remaining_nodes[i];
        assert( 0 <= current_node and current_node < counts[dim] );
        
        // What is the dimension of the subsimplex around which the interface is made?        
        int k = dim - how_many_connected_faces[i];
        Assert( i == 0 or k <= dim, i, k );

        if( k == dim ) assert( current_prefix.size() == 0 );
        if( current_prefix.size() == 0 ) assert( k == dim );

        
        // Determine whether the new simplex can be part of a shelling,
        // and if yes, what is the common subsimplex

        int common_subsimplex = mesh.nullindex; 

        {

            // get the subsimplices of the proposed node of dimension k
            const auto k_subsimplices_of_node = mesh.getsubsimplices( dim, k, current_node ).getvalues();

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

            // check that all parents of the k subsimplex are here
            const auto& parents = mesh.getsupersimplices( dim, k, common_subsimplex );
            std::vector<bool> parent_is_already_here( parents.size(), false );

            for( int p = 0; p < parents.size(); p++ )
            {
                if( parents[p] == current_node ) parent_is_already_here[p] = true;

                for( int s : current_prefix ) { 
                    if( parents[p] == s ) parent_is_already_here[p] = true;
                }
            }

            bool all_parents_are_here = std::all_of( parent_is_already_here.begin(), parent_is_already_here.end(), [](bool b){return b;} );

            shelling_compatible[i] = all_parents_are_here;
        }



        // Having determined whether the next volume is compatible, 
        // let us now compute the next vector in the coefficient table
        {

            Float PF = 1. / Constants::pi ; 

            Float pullbackfactor_d = ( k != dim ? info.C5[current_node][k][form_degree+1] : 0. );
            Float pullbackfactor_0 = ( k != dim ? info.C5[current_node][k][form_degree  ] : 0. );
            
            Float A = PF;
            Float B = PF * pullbackfactor_d;
            Float C =      pullbackfactor_0;

            // obtain all previous indices of the common subsimplex of dimension n-k
            const auto& relevant_volumes = mesh.getsupersimplices( dim, k, common_subsimplex );

            FloatVector new_coefficients( counts[dim], 0. );

            new_coefficients[current_node] = A;

            // for( auto p : relevant_volumes ) LOG << p << space; LOG << nl;

            for( int j = 0; j < current_prefix.size(); j++ )
            {
                int previous_node = current_prefix[j];

                bool is_relevant = std::find( relevant_volumes.begin(), relevant_volumes.end(), previous_node ) != relevant_volumes.end();

                if( not is_relevant ) continue;

                new_coefficients[j] += B;

                new_coefficients += C * coefficient_table[j];
            }

            coefficient_vectors[i] = new_coefficients;
            
            weight_for_node_reflect[i] = new_coefficients.l2norm();

            weight_for_node_morphin[i] = ( k != dim ? info.C7[current_node][k][form_degree] * info.C8[current_node][k][form_degree+1] : PF ); // TODO: which form degree?

        }
    
    }

    
    
    // unless we are at the start:
    // delete all the nodes that don't give a shelling 
    if( current_prefix.size() > 0 )
    {
        int i = 0;
        while( i < remaining_nodes.size() )
        {
            if( shelling_compatible[i] ) {
                assert( not std::isnan( weight_for_node_reflect[i] ) );
                i++;
            } else {
                remaining_nodes.erase(          remaining_nodes.begin() + i          );
                face_is_connected.erase(        face_is_connected.begin() + i        );
                how_many_connected_faces.erase( how_many_connected_faces.begin() + i );
                shared_subsimplex_dim.erase(    shared_subsimplex_dim.begin() + i    );
                shelling_compatible.erase(      shelling_compatible.begin() + i      );
                coefficient_vectors.erase(      coefficient_vectors.begin() + i      );
                weight_for_node_reflect.erase(  weight_for_node_reflect.begin() + i  );
                weight_for_node_morphin.erase(  weight_for_node_morphin.begin() + i  );
            }
        }

        assert( remaining_nodes.size() > 0 );
        assert( remaining_nodes.size() == face_is_connected.size()        );
        assert( remaining_nodes.size() == how_many_connected_faces.size() );
        assert( remaining_nodes.size() == shared_subsimplex_dim.size()    );
        assert( remaining_nodes.size() == shelling_compatible.size()      );
        assert( remaining_nodes.size() == coefficient_vectors.size()      );
        assert( remaining_nodes.size() == weight_for_node_reflect.size()  );
        assert( remaining_nodes.size() == weight_for_node_morphin.size()  );
    }

    // unless we are at the start:
    // sort nodes by priority
    // simply some lazy bubble sort
    /////if( current_prefix.size() > 0 )
    for( int i = 0; i < remaining_nodes.size(); i++ )
    for( int j = 1; j < remaining_nodes.size(); j++ )
    {
        bool correct_order = weight_for_node_reflect[j-1] <= weight_for_node_reflect[j];

        if( correct_order ) continue;

        std::swap( remaining_nodes[j-1],          remaining_nodes[j]          );
        std::swap( how_many_connected_faces[j-1], how_many_connected_faces[j] );
        std::swap( face_is_connected[j-1],        face_is_connected[j]        );
        std::swap( shared_subsimplex_dim[j-1],    shared_subsimplex_dim[j]    );
        //std::swap( shelling_compatible[j-1],      shelling_compatible[j]      );
        {
            bool temp = shelling_compatible[j-1];
            shelling_compatible[j-1] = shelling_compatible[i];
            shelling_compatible[i] = temp;
        }
        std::swap(     coefficient_vectors[j-1],     coefficient_vectors[j] );
        std::swap( weight_for_node_reflect[j-1], weight_for_node_reflect[j] );
        std::swap( weight_for_node_morphin[j-1], weight_for_node_morphin[j] );

        Assert( weight_for_node_reflect[j-1] <= weight_for_node_reflect[j], weight_for_node_reflect[j-1], weight_for_node_reflect[j] );
    }

    for( int i = 1; i < remaining_nodes.size(); i++ ) Assert( weight_for_node_reflect[i] >= weight_for_node_reflect[i-1], weight_for_node_reflect[i], weight_for_node_reflect[i-1] );
    

    // only reachable nodes are left and they are ordered by priority 

    for( int i = 0; i < remaining_nodes.size(); i++ )
    {
        int node = remaining_nodes[i];
        
        auto next_prefix = current_prefix;
        next_prefix.push_back( node );

        auto next_coefficient_table = coefficient_table;
        next_coefficient_table.push_back( coefficient_vectors[i] );

        assert( next_prefix.size() == next_coefficient_table.size() );

        Float next_multweight = multweight * weight_for_node_morphin[i];

        generate_shellings2( mesh, form_degree, info, shellings_found, next_prefix, next_coefficient_table, next_multweight );
    }
    
}













#endif // SHELLING_LISTER_2_HPP
