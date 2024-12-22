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

    std::vector<std::vector<FloatVector>> midpoints;                // [dimension][simplex index]
    std::vector<std::vector<Float>>       heights_of_vertex;        // [volume index][local vertex index]
    std::vector<std::vector<FloatVector>> heightvectors_of_vertex;  // [volume index][local vertex index]
    
    std::vector<std::vector<std::vector<std::vector<Float>>>> C5;
    std::vector<std::vector<std::vector<std::vector<Float>>>> C6;
    std::vector<std::vector<std::vector<std::vector<Float>>>> C7;
    std::vector<std::vector<std::vector<std::vector<Float>>>> C8;
    
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

        const auto faces = mesh.get_subsimplices(dim,dim-1,i).getvalues(); assert( faces.size() == dim+1 );

        std::vector<Float> heights( dim+1 );
        for( int face_index = 0; face_index < faces.size(); face_index++ ) heights[face_index] = dim * volumes[i] / mesh.getMeasure( dim-1, faces[face_index] );

        Float min_height = *min_element( heights.begin(), heights.end() );

        min_height_of_vertex = minimum( min_height, min_height_of_vertex );

        for( int face_index = 0; face_index <= dim; face_index++ ) {
            int vi = mesh.get_opposite_subsimplex_index( dim, dim-1, i, face_index );
            heights_of_vertex[i][vi] = heights[face_index];
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
    for( int vertex_index = 0; vertex_index <= dim; vertex_index++ )
    {
        heightvectors_of_vertex[i][vertex_index] = mesh.getHeightVector( dim, i, vertex_index );
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
    for( int face = 0; face < counts[dim-1]; face++ )
    {
        const auto parents = mesh.get_supersimplices(dim,dim-1,face);
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

    C5.resize( counts[dim], std::vector<std::vector<std::vector<Float>>>( dim+1, std::vector<std::vector<Float>>( dim, std::vector<Float>() ) ) );
    C6.resize( counts[dim], std::vector<std::vector<std::vector<Float>>>( dim+1, std::vector<std::vector<Float>>( dim, std::vector<Float>() ) ) );
    C7.resize( counts[dim], std::vector<std::vector<std::vector<Float>>>( dim+1, std::vector<std::vector<Float>>( dim, std::vector<Float>() ) ) );
    C8.resize( counts[dim], std::vector<std::vector<std::vector<Float>>>( dim+1, std::vector<std::vector<Float>>( dim, std::vector<Float>() ) ) );

    // Compute the coefficients 
    for( int i = 0; i < counts[dim]; i++ ) 
    for( int form_degree = 0; form_degree <= dim; form_degree++ )
    for( int level = 0; level < dim; level++ )
    {
        const int num_of_subsimplices = binomial_integer(dim+1,level+1);

        C5[i][form_degree][level] = std::vector<Float>( num_of_subsimplices );
        C6[i][form_degree][level] = std::vector<Float>( num_of_subsimplices );
        C7[i][form_degree][level] = std::vector<Float>( num_of_subsimplices );
        C8[i][form_degree][level] = std::vector<Float>( num_of_subsimplices );
        
        for( int sub_index = 0; sub_index < num_of_subsimplices; sub_index++ )
        {
        
            Float singular_max     = 0.;
            Float singular_min_inv = 0.;
            Float singular_prod    = 0.;

            Float tmax    = 0.;
            Float tinvmin = 0.;
            Float tdet    = 0.;
            
            const int n = dim;

            // rough initial estimates 
            {
                const Float Ctheta = max_diameter_ratio;

                const Float kappa  = aspect_condition_number[i];

                const Float rho    = minimum( 1., min_height_of_vertex / (level+1.) );

                const Float mu_l   = 0.5 * sqrt( square(1+rho) + square( square(1+rho) * (n-level) * kappa ) ) + 0.5 * sqrt( square(1-rho) + square( square(1+rho) * (n-level) * kappa ) );
                
                const Float xi_l   = sqrt( square(2*rho+1) + square( (n-level) * kappa ) ) + sqrt( 1. + square( (n-level) * kappa ) );

                assert( std::isfinite( Ctheta ) && Ctheta > 0 );
                assert( std::isfinite( kappa  ) && kappa  > 0 );
                assert( std::isfinite( mu_l   ) && mu_l   > 0 );

                singular_max     = mu_l;
                singular_min_inv = mu_l / rho;
                singular_prod    = rho;
            
                tmax    = xi_l / ( 2 * (1.+rho) );
                tinvmin = xi_l / ( 2 * rho      );
                tdet    = rho / ( 1 + rho );
                
            }

            // improved estimates, simplex by simplex 
            // if(false)
            {
                // find the radius that we can use

                Float rho = 1.;

                const int k = level;

                int subsimplex = mesh.get_subsimplex( dim, k, i, sub_index );

                auto parents  = mesh.get_supersimplices( dim, k, subsimplex );
                auto vertices = mesh.get_subsimplices  (   k, 0, subsimplex ).getvalues();

                for( int vertex : vertices )
                for( int parent : parents  )
                {
                    int vertex_localindex = mesh.get_subsimplex_index( dim, 0, parent, vertex );
                    Float height = mesh.getHeight( dim, parent, vertex_localindex );
                    rho = minimum( rho, height / (k+1) );
                }

                // find the heights 

                int index_oppsite_subsimplex = mesh.get_opposite_subsimplex_index( dim, k, i, sub_index );
                int opposite_subsimplex      = mesh.get_subsimplex( dim, dim - k - 1, i, index_oppsite_subsimplex );

                auto other_vertices = mesh.get_subsimplices( dim-k-1, 0, opposite_subsimplex ).getvalues();
                
                auto zs = mesh.get_midpoint(       k,          subsimplex );
                auto zc = mesh.get_midpoint( dim-k-1, opposite_subsimplex );

                std::vector<FloatVector> height_vectors( dim-k, FloatVector(dim,notanumber) );

                for( int vi = 0; vi <= dim-k-1; vi++ ) {
                    height_vectors[vi] = mesh.getHeightVector( dim, i, vi );
                    LOG << height_vectors[vi].scalarproductwith( zc - zs ) << tab;
                } LOG << nl;

                // computation 
                // when S has got k+1 vertices, then S' must have dim-k vertices 
                // Then S has got dimension k and S' has got dimension dim-k-1
                
                Float improved_singular_max = 0.; 

                Float matrixwise_singular_max = 0.; 
                
                for( int vi = 0; vi <= dim-k-1; vi++ ) 
                {
                    
                    auto z = zc - zs;
                    auto h = height_vectors[vi] / ( dim-k );
                    auto b = z - h; //( z.scalarproductwith(h) / h.norm_sq() ) * h;
                    Assert( ( z - h - b ).is_numerically_small(), z, '\n', h, '\n', b );
                    auto tanbeta = b.norm() / h.norm();

                    Float proposed_improved_singular_max 
                    = 
                    0.5 * sqrt( square( 1. + rho ) + square( ( 1. + rho ) * tanbeta ) ) 
                    +
                    0.5 * sqrt( square( 1. - rho ) + square( ( 1. + rho ) * tanbeta ) );

                    improved_singular_max = maximum( improved_singular_max, proposed_improved_singular_max );
                }

                for( int v1 = 0; v1 <=       k; v1++ ) 
                for( int v2 = 0; v2 <= dim-k-1; v2++ ) 
                {

                    std::vector<FloatVector> columns;
                    for( int w1 = 0; w1 <= k; w1++ ) 
                        if( v1 != w1 ) 
                            columns.push_back( mesh.getCoordinates().getvectorclone(       vertices[w1] ) - zs );

                    for( int w2 = 0; w2 <= dim-k-1; w2++ ) 
                        if( v2 != w2 ) 
                            columns.push_back( mesh.getCoordinates().getvectorclone( other_vertices[w2] ) - zs );

                    columns.push_back( zc - zs );

                    assert( columns.size() == dim );

                    DenseMatrix trafo1( dim, dim, columns );

                    auto trafo1inv = Inverse( trafo1 );

                    columns.back().scale( -rho );

                    DenseMatrix trafo2( dim, dim, columns );

                    auto fulltrafo = trafo2 * trafo1inv;

                    Float proposed_matrixwise_singular_max = fulltrafo.operator_norm_estimate();

                    matrixwise_ ( matrixwise_singular_max, proposed_matrixwise_singular_max );

                    // LOG << "DET " << Determinant( fulltrafo ) / (-rho) << nl;
                    // LOG << "MAX " << fulltrafo.operator_norm_estimate() / singular_max << nl;
                    
                }
                
                assert( improved_singular_max > 0. );
                
                assert( matrixwise_singular_max > 0. );

                LOG << "SMAX " << improved_singular_max / matrixwise_singular_max << nl;

                singular_max = minimum( singular_max, improved_singular_max, matrixwise_singular_max );
                
                assert( singular_max > 0. ); 

                singular_prod = rho;
                singular_min_inv = singular_max / rho;

                PING;
                

            }


            // slight improvement in the case of faces 
            // if(false)
            if( level == dim-1 ) {
                
                const Float Ctheta = max_diameter_ratio;

                const Float kappa  = aspect_condition_number[i];

                Float facebased_singular_max = 0.5 * sqrt( square( Ctheta * kappa + 1. ) + kappa ) + 0.5 * sqrt( square( Ctheta * kappa - 1. ) + kappa );

                singular_max = minimum( singular_max, facebased_singular_max );

                singular_min_inv = singular_max / max_volume_ratio;
                singular_prod = max_volume_ratio;
            }
            
            // better improvement in the case of faces 
            // if(false)
            if( level == dim-1 ) {
                
                assert( 0 <= sub_index and sub_index < binomial_integer(dim+1,level+1) );
                assert( 0 <= sub_index and sub_index < binomial_integer(dim+1,dim) );
                    
                int face = mesh.get_subsimplex( dim, dim-1, i, sub_index ); // TODO: Debug 

                assert( 0 <= face and face < counts[dim-1] );

                auto parents = mesh.get_supersimplices(dim,dim-1,face);
                
                if( parents.size() == 2 ) 
                {
                    assert( parents[0] == i or parents[1] == i );

                    int other_parent = ( parents[0] == i ? parents[1] : parents[0] );
                    
                    assert( other_parent != i and 0 <= other_parent and other_parent < counts[dim-1] );

                    int localindex_face_first = sub_index;
                    int localindex_face_other = mesh.get_subsimplex_index(dim,dim-1,other_parent,face);

                    int opposite_vertex_localindex_first = mesh.get_opposite_subsimplex_index( dim, dim-1,            i, localindex_face_first );
                    int opposite_vertex_localindex_other = mesh.get_opposite_subsimplex_index( dim, dim-1, other_parent, localindex_face_other );

                    int vertex_first = mesh.get_subsimplex( dim, 0,            i, opposite_vertex_localindex_first );
                    int vertex_other = mesh.get_subsimplex( dim, 0, other_parent, opposite_vertex_localindex_other );

                    auto height_vector_first = mesh.getHeightVector( dim,            i, opposite_vertex_localindex_first );
                    auto height_vector_other = mesh.getHeightVector( dim, other_parent, opposite_vertex_localindex_other );

                    auto midpoint_of_face = mesh.get_midpoint( dim-1, face );

                    Float a = height_vector_other.norm() /  height_vector_first.norm();

                    LOG << vertex_first << space << vertex_other << nl;

                    auto z_first = mesh.getCoordinates().getvectorclone( vertex_first );
                    auto z_other = mesh.getCoordinates().getvectorclone( vertex_other );

                    auto b_first = z_first - z_first.scalarproductwith(height_vector_first) / height_vector_first.norm_sq() * height_vector_first;
                    auto b_other = z_other - z_other.scalarproductwith(height_vector_other) / height_vector_other.norm_sq() * height_vector_other;

                    Float c = ( b_first - b_other ).norm() / height_vector_first.norm();

                    Float facebased_singular_max = 0.5 * sqrt( square( a + 1. ) + c ) + 0.5 * sqrt( square( a - 1. ) + c );
                    

                    singular_max = minimum( singular_max, facebased_singular_max );

                    singular_min_inv = singular_max / a;
                    singular_prod = a;

                } else {

                    assert( parents.size() == 1 );
                    assert( parents[0] == i );

                }

            }
            
            if( form_degree == 0 ) {
                C5[i][form_degree][level][sub_index] = power_numerical( singular_prod, -n/2. );
                C6[i][form_degree][level][sub_index] = power_numerical( singular_prod, +n/2. );
            } else if( 0 < form_degree && form_degree < dim ) {
                C5[i][form_degree][level][sub_index] = singular_max    * power_numerical( singular_prod, -n/2. );
                C6[i][form_degree][level][sub_index] = singular_min_inv * power_numerical( singular_prod, +n/2. );
            } else {
                assert( form_degree == dim );
                C5[i][form_degree][level][sub_index] = power_numerical( singular_prod, +1.-n/2. );
                C6[i][form_degree][level][sub_index] = power_numerical( singular_prod, -1.+n/2. );
            }

            if( form_degree == 0 ) {
                C7[i][form_degree][level][sub_index] = power_numerical( tdet, -n/2. );
                C8[i][form_degree][level][sub_index] = power_numerical( tdet, +n/2. );
            } else if( 0 < form_degree && form_degree < dim ) {
                C7[i][form_degree][level][sub_index] = tmax    * power_numerical( tdet, -n/2. );
                C8[i][form_degree][level][sub_index] = tinvmin * power_numerical( tdet, +n/2. );
            } else {
                assert( form_degree == dim );
                C7[i][form_degree][level][sub_index] = power_numerical( tdet, +1.-n/2. );
                C8[i][form_degree][level][sub_index] = power_numerical( tdet, -1.+n/2. );
            }

            Assert( not std::isnan( C5[i][form_degree][level][sub_index] ), form_degree ); Assert( C5[i][form_degree][level][sub_index] > 0. ); 
            Assert( not std::isnan( C6[i][form_degree][level][sub_index] ), form_degree ); Assert( C6[i][form_degree][level][sub_index] > 0. ); 
            Assert( not std::isnan( C7[i][form_degree][level][sub_index] ), form_degree ); Assert( C7[i][form_degree][level][sub_index] > 0. ); 
            Assert( not std::isnan( C8[i][form_degree][level][sub_index] ), form_degree ); Assert( C8[i][form_degree][level][sub_index] > 0. ); 
            
        }

    }

    for( auto a : C5 ) for( auto b : a ) for( int i = 0; i < dim; i++ ) assert( b[i].size() == binomial_integer(dim+1,i+1) );

    for( auto a : C5 ) for( auto b : a ) for( auto c : b ) for( auto d : c ) assert( not std::isnan(d) && d > 0. );
    for( auto a : C6 ) for( auto b : a ) for( auto c : b ) for( auto d : c ) assert( not std::isnan(d) && d > 0. );
    for( auto a : C7 ) for( auto b : a ) for( auto c : b ) for( auto d : c ) assert( not std::isnan(d) && d > 0. );
    for( auto a : C8 ) for( auto b : a ) for( auto c : b ) for( auto d : c ) assert( not std::isnan(d) && d > 0. );
    
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

    const bool sort_nodes_by_priority = false;

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

        LOG << "FOUND: ";
        
        Float estimate1 = 0.;
        for( int i = 0; i < counts[dim]; i++ ) {
            estimate1 += coefficient_table[i].norm_sq();
            LOGPRINTF( "%i %.4f\t", current_prefix[i], estimate1 );
        }
        estimate1 = sqrt(estimate1);

        LOGPRINTF( "w1=%e w2=%e\n", estimate1, multweight );

        auto insertion_point = shellings_found.begin();
        while( insertion_point != shellings_found.end() && insertion_point->weight_reflection < estimate1 ) insertion_point++;
        shellings_found.insert( insertion_point, shelling( current_prefix, estimate1, multweight ) );

        if(false)
        for( int i = 0; i < shellings_found.size(); i++ ) {
            //assert( shellings_found[i-1].weight_reflection <= shellings_found[i].weight_reflection );
            LOG << tab << shellings_found[i].weight_reflection << nl;
        }
            
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
        
        const auto faces_of_node = mesh.get_subsimplices( dim, dim-1, node ).getvalues();

        assert( faces_of_node.size() == face_is_connected[i].size() );

        how_many_connected_faces[i] = 0;

        for( int index_f = 0; index_f < faces_of_node.size(); index_f++ )
        {
            const int face = faces_of_node[index_f];
            
            const auto parents = mesh.get_supersimplices( dim, dim-1, face );
            
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


    


    
    
    
    // compute the weight associated with the extension along each possible node 

    assert( remaining_nodes.size() > 0 );

    std::vector<int>   shared_subsimplex_dim( remaining_nodes.size(), -1    );
    std::vector<int>   shared_subsimplex    ( remaining_nodes.size(), mesh.nullindex );

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

        shared_subsimplex_dim[i] = k;

        if( k == dim ) assert( current_prefix.size() == 0 );
        if( current_prefix.size() == 0 ) assert( k == dim );

        
        // Determine whether the new simplex can be part of a shelling,
        // and if yes, what is the common subsimplex

        int common_subsimplex = mesh.nullindex; 

        {

            // get the subsimplices of the proposed node of dimension k
            const auto k_subsimplices_of_node = mesh.get_subsimplices( dim, k, current_node ).getvalues();

            const auto faces_of_node = mesh.get_subsimplices( dim, dim-1, current_node ).getvalues();

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

            {
                assert( 0 <= common_subsimplex and common_subsimplex < mesh.count_simplices(k) );
                const int index = mesh.get_subsimplex_index( dim, k, current_node, common_subsimplex );
                assert( 0 <= index and index < mesh.count_subsimplices(dim,k) );
                assert( mesh.get_subsimplices(dim,k,current_node)[index] == common_subsimplex );
            }

            // check that all parents of the k subsimplex are here
            const auto& parents = mesh.get_supersimplices( dim, k, common_subsimplex );
            std::vector<bool> parent_is_already_here( parents.size(), false );

            for( int p = 0; p < parents.size(); p++ )
            {
                if( parents[p] == current_node ) parent_is_already_here[p] = true;

                for( int s : current_prefix ) { 
                    if( parents[p] == s ) parent_is_already_here[p] = true;
                }
            }

            bool all_parents_are_here = std::all_of( parent_is_already_here.begin(), parent_is_already_here.end(), [](bool b){return b;} );

            shared_subsimplex[i] = ( all_parents_are_here ? common_subsimplex : mesh.nullindex );
        }



        // Having determined whether the next volume is compatible, 
        // let us now compute the next vector in the coefficient table
        if( shared_subsimplex[i] != mesh.nullindex )
        {
            assert( 0 <= shared_subsimplex[i] and shared_subsimplex[i] < counts[k] );

            int sub_index = 0;
            // const int binom = binomial_integer(dim+1,k+1);
            // while( sub_index < binom ) { if( mesh.get_subsimplices(dim,k,current_node)[sub_index] == shared_subsimplex[current_node] ) break; sub_index++; } assert( sub_index < binom ); 
            sub_index = mesh.get_subsimplex_index( dim, k, current_node, shared_subsimplex[i] );
            assert( mesh.get_subsimplices(dim,k,current_node)[sub_index] == shared_subsimplex[i] );

            Float PF = 2. / Constants::pi ; // needs to be made rigorous

            Float pullbackfactor_d = ( k != dim ? info.C5[current_node][form_degree+1][k][sub_index] : 0. );
            Float pullbackfactor_0 = ( k != dim ? info.C5[current_node][form_degree  ][k][sub_index] : 0. );
            
            Float A = PF;
            Float B = PF * pullbackfactor_d;
            Float C =      pullbackfactor_0;

            // obtain all previous indices of the common subsimplex of dimension n-k
            const auto& relevant_volumes = mesh.get_supersimplices( dim, k, common_subsimplex );

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

            weight_for_node_morphin[i] = ( k != dim ? info.C7[current_node][form_degree][k][sub_index] * info.C8[current_node][form_degree+1][k][sub_index] : PF ); // TODO: which form degree?

        }
    
    }

    
    
    // unless we are at the start:
    // delete all the nodes that don't give a shelling 
    if( current_prefix.size() > 0 )
    {
        int i = 0;
        while( i < remaining_nodes.size() )
        {
            if( shared_subsimplex[i] != mesh.nullindex ) {
                assert( how_many_connected_faces[i] > 0 );
                assert( dim - how_many_connected_faces[i] == shared_subsimplex_dim[i] );
                assert( not std::isnan( weight_for_node_reflect[i] ) );
                i++;
            } else {
                remaining_nodes.erase(          remaining_nodes.begin() + i          );
                face_is_connected.erase(        face_is_connected.begin() + i        );
                how_many_connected_faces.erase( how_many_connected_faces.begin() + i );
                shared_subsimplex_dim.erase(    shared_subsimplex_dim.begin() + i    );
                shared_subsimplex.erase(        shared_subsimplex.begin() + i        );
                coefficient_vectors.erase(      coefficient_vectors.begin() + i      );
                weight_for_node_reflect.erase(  weight_for_node_reflect.begin() + i  );
                weight_for_node_morphin.erase(  weight_for_node_morphin.begin() + i  );
            }
        }

        assert( remaining_nodes.size() == face_is_connected.size()        );
        assert( remaining_nodes.size() == how_many_connected_faces.size() );
        assert( remaining_nodes.size() == shared_subsimplex_dim.size()    );
        assert( remaining_nodes.size() == shared_subsimplex.size()        );
        assert( remaining_nodes.size() == coefficient_vectors.size()      );
        assert( remaining_nodes.size() == weight_for_node_reflect.size()  );
        assert( remaining_nodes.size() == weight_for_node_morphin.size()  );

        if( current_prefix.size () == 0 ) assert( remaining_nodes.size() > 0 );        
    }

    for( int j = 0; j < remaining_nodes.size(); j++ )
    {
        assert( coefficient_vectors[j].isfinite()           );
        assert( std::isfinite( weight_for_node_reflect[j] ) );
        assert( std::isfinite( weight_for_node_morphin[j] ) );
    }


    if( remaining_nodes.size() == 0 ) return;

    if(false)
    {
        LOG << "current nodes:   "; for( auto node : current_prefix  ) LOG << node << space; 
        LOG << tab << tab << shellings_found.size();
        if( shellings_found.size() > 0 ) LOG << tab << shellings_found.front().weight_reflection << "--" << shellings_found.back().weight_reflection << tab;
        LOG << nl;
        LOG << "remaining nodes: "; for( auto node : remaining_nodes ) LOG << node << space; LOG << nl;
    }




    // unless we are at the start:
    // sort nodes by priority
    // simply some lazy bubble sort
    /////if( current_prefix.size() > 0 )
    if( sort_nodes_by_priority )
    for( int i = 0; i < remaining_nodes.size(); i++ )
    for( int j = 1; j < remaining_nodes.size(); j++ )
    {
        bool correct_order = ( false and how_many_connected_faces[j-1] < how_many_connected_faces[j] ) or ( weight_for_node_reflect[j-1] <= weight_for_node_reflect[j] );

        if( correct_order ) continue;

        std::swap( remaining_nodes[j-1],          remaining_nodes[j]          );
        std::swap( how_many_connected_faces[j-1], how_many_connected_faces[j] );
        std::swap( face_is_connected[j-1],        face_is_connected[j]        );
        std::swap( shared_subsimplex_dim[j-1],    shared_subsimplex_dim[j]    );
        std::swap( shared_subsimplex[j-1],        shared_subsimplex[j]        );
        // { bool temp = shared_subsimplex[j-1]; shared_subsimplex[j-1] = shared_subsimplex[i]; shared_subsimplex[i] = temp; }
        std::swap(     coefficient_vectors[j-1],     coefficient_vectors[j] );
        std::swap( weight_for_node_reflect[j-1], weight_for_node_reflect[j] );
        std::swap( weight_for_node_morphin[j-1], weight_for_node_morphin[j] );

        Assert( ( false and how_many_connected_faces[j-1] < how_many_connected_faces[j] ) or ( weight_for_node_reflect[j-1] <= weight_for_node_reflect[j] ) );
    }

    if( sort_nodes_by_priority )
    for( int i = 1; i < remaining_nodes.size(); i++ ) 
        Assert( ( false and how_many_connected_faces[i-1] < how_many_connected_faces[i] ) or weight_for_node_reflect[i-1] <= weight_for_node_reflect[i], weight_for_node_reflect[i], weight_for_node_reflect[i-1] );
    
    // only reachable nodes are left and they are ordered by priority 

    Float weight_reflection_so_far = 0.;
    for( auto coefficient : coefficient_table ) weight_reflection_so_far += coefficient.norm_sq();

    for( int i = 0; i < remaining_nodes.size(); i++ )
    {
        Float weight_here = coefficient_vectors[i].norm_sq();

        if( shellings_found.size() > 0 and sqrt( weight_here + weight_reflection_so_far ) >= shellings_found.front().weight_reflection ) 
            if( sort_nodes_by_priority )
                break;
            else 
                continue;
        
        
        int node = remaining_nodes[i];
        
        auto next_prefix = current_prefix;
        next_prefix.push_back( node );

        auto next_coefficient_table = coefficient_table;
        next_coefficient_table.push_back( coefficient_vectors[i] );

        assert( next_prefix.size() == next_coefficient_table.size() );

        Float next_multweight = multweight * weight_for_node_morphin[i];

        generate_shellings2( mesh, form_degree, info, shellings_found, next_prefix, next_coefficient_table, next_multweight );
    }

    static int call_counter = 0;
    call_counter++;
    if( current_prefix.size() == 0 ) LOG << "Calls: " << call_counter << nl;
    
}













#endif // SHELLING_LISTER_2_HPP
