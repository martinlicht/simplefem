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
    
    std::vector<std::vector<Float>> C5;
    std::vector<std::vector<Float>> C5prime;
    std::vector<std::vector<Float>> C5det;
    std::vector<std::vector<Float>> C6;
    std::vector<std::vector<Float>> C6prime;
    std::vector<std::vector<Float>> C6det;
    std::vector<std::vector<Float>> C7;
    std::vector<std::vector<Float>> C7prime;
    std::vector<std::vector<Float>> C7det;
    std::vector<std::vector<Float>> C8;
    std::vector<std::vector<Float>> C8prime;
    std::vector<std::vector<Float>> C8det;

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
        
        C5.resize(counts[dim]+1,      std::vector<Float>(counts[dim],notanumber) );
        C5prime.resize(counts[dim]+1, std::vector<Float>(counts[dim],notanumber) );
        C5det.resize(counts[dim]+1,   std::vector<Float>(counts[dim],notanumber) );
        C6.resize(counts[dim]+1,      std::vector<Float>(counts[dim],notanumber) );
        C6prime.resize(counts[dim]+1, std::vector<Float>(counts[dim],notanumber) );
        C6det.resize(counts[dim]+1,   std::vector<Float>(counts[dim],notanumber) );
        C7.resize(counts[dim]+1,      std::vector<Float>(counts[dim],notanumber) );
        C7prime.resize(counts[dim]+1, std::vector<Float>(counts[dim],notanumber) );
        C7det.resize(counts[dim]+1,   std::vector<Float>(counts[dim],notanumber) );
        C8.resize(counts[dim]+1,      std::vector<Float>(counts[dim],notanumber) );
        C8prime.resize(counts[dim]+1, std::vector<Float>(counts[dim],notanumber) );
        C8det.resize(counts[dim]+1,   std::vector<Float>(counts[dim],notanumber) );
        
        // estimate the algebraic condition number of these guys
    
        for( int i = 0; i < counts[dim]; i++ ) 
        {
            diameters[i] = mesh.getDiameter(dim,i);
            
            volumes[i] = mesh.getMeasure(dim,i);

            const auto faces = mesh.getsubsimplices(dim,dim-1,i).getvalues();

            std::vector<Float> heights( dim+1 );
            for( int f = 0; f < faces.size(); f++ ) 
                heights[f] = dim * volumes[i] / mesh.getMeasure( dim-1, faces[f] );

            Float height = *min_element( heights.begin(), heights.end() );

            aspect_condition_number[i] = diameters[i] / height;

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
        for( int k = 0; k <= dim; k++ )
        {
            const   int n = dim;

            const Float B = sqrt( 1. + square( (n-k)*aspect_condition_number[i] ) ) - (n-k)*aspect_condition_number[i];

            const Float Ctheta = max_diameter_ratio;

            const Float kappa = aspect_condition_number[i];

            LOG << i << ":\t" << k << ":\t" << B << space << Ctheta << space << kappa << nl;


            const Float Psi_estimate = sqrt( 1. + square( (k+1) * (1+B) * kappa ) / 4 ) + ( (k+1) * (1+B) * kappa ) / 2;
            
            C5[k][i]      = (k+1) * Ctheta * kappa * (B*B) * Psi_estimate;
            
            C5prime[k][i] = (k+1) * B * Ctheta * kappa * B;

            C5det[k][i]   = power_numerical( (k+1) * B * Ctheta * kappa, n );

            
            C6[k][i]      = (k+1) * Ctheta * kappa * B * Psi_estimate;
            
            C6prime[k][i] = (k+1) * B * Ctheta * kappa;

            C6det[k][i]   = power_numerical( (k+1) * Ctheta * kappa, n );

            
            C7[k][i]      = ( 1. + 3./2. * (k+1) * kappa ) * (k+1) * kappa * B * Psi_estimate;
            
            C7prime[k][i] = ( 1. + 3./2. * (k+1) * kappa ) * (k+1) * B * Ctheta * kappa;

            C7det[k][i]   = 0.5 * power_numerical( (k+1) * B * Ctheta * kappa, n );


            C8[k][i]      = ( 2. + 3. * (k+1) * kappa ) * (k+1) * kappa * Psi_estimate;
            
            C8prime[k][i] = ( 2. + 3. * (k+1) * kappa ) * (k+1) * Ctheta * kappa;

            C8det[k][i]   = 2.0 * power_numerical( (k+1) * Ctheta * kappa, n );
        }

        

    }
};








Float estimate_shelling_quality( 
    const Mesh& mesh,
    std::vector<int> shelling,
){
    
    

    

    


    std::vector<FloatVector> coefficient_table( counts[dim], FloatVector( counts[dim], 0. ) );
    
    // start computing the relevant data on each piece 
    for( int i = 0; i < counts[dim]; i++ )
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
