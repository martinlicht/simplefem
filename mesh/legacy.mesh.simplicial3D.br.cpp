
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"

#include "mesh.simplicial3D.hpp"















/*
 * * * * * VERTEX BISECTION
 */
    
void MeshSimplicial3D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    check();
    
    /*******************/
    /* DATA COLLECTION */
    /*******************/
    
    // gather basic data about the bisection, in particular the number of simplices involved.
    
    int e_back_vertex  = data_edge_vertices[ e ][ 0 ];
    int e_front_vertex = data_edge_vertices[ e ][ 1 ];
    
    std::vector<int> old_faces      = get_face_parents_of_edge( e );
    std::vector<int> old_tetrahedra = get_tetrahedron_parents_of_edge( e );
    
    std::vector<int> localindex_of_face_refinementedge( old_faces.size() );
    std::vector<int> localindex_of_tetrahedron_refinementedge( old_tetrahedra.size() );
    
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    
    
    /********************************/
    /* ALLOCATE MEMORY FOR THE DATA */
    /********************************/
    
    // allocate additional memory for the new simplex data after bisection 
    
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_faces.resize               ( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_face_firstparent_tetrahedron.resize    ( counter_faces + old_faces.size() + old_tetrahedra.size(),                                      nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_edges.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_edge_firstparent_tetrahedron.resize    ( counter_edges + old_faces.size() + 1,                                                                  nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_vertices.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                                            nullindex );
    
    data_face_nextparents_of_edges.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_edges.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_edge_firstparent_face.resize    ( counter_edges     + old_faces.size() + 1,                                           nullindex );
    
    data_face_nextparents_of_vertices.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_vertices.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_face.resize     ( counter_vertices + 1,                                                               nullindex );
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_edge_vertices.resize               ( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,                                   nullindex );
    
    
    
    
    /****************/
    /* FILL IN DATA */
    /****************/
    
    /* vertices of the children of the bisected edge */
    
    data_edge_vertices[ e             ][ 0 ] = e_back_vertex;
    data_edge_vertices[ e             ][ 1 ] = counter_vertices;
      
    data_edge_vertices[ counter_edges ][ 0 ] = counter_vertices;
    data_edge_vertices[ counter_edges ][ 1 ] = e_front_vertex;
    
    /* next parent of back vertex in back edge is already set */
    /* next parent of front vertex in front edge is inherited */
    
    data_edge_nextparents_of_vertices[ counter_edges ][ 1 ] = data_edge_nextparents_of_vertices[ e ][ 1 ];
    
    // edge parent list of back vertex stays the same 
    // run over the front vertex edge parent list and replace 'e' by 'counter_edges'
    
    if( data_vertex_firstparent_edge[ e_front_vertex ] == e ) {
      
      // if 'e' is the first edge parent of the front vertex, just replace that by 'counter_edges'
      data_vertex_firstparent_edge[ e_front_vertex ] = counter_edges;
      
    } else {
      
      // else, find 'e' in the edge parent list of the front vertex and replace 
      int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
      
      while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
        current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
      
      data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edges;
      
    }
    
    /* first and second parent of new vertex are the children of the bisected edge */
    
    data_vertex_firstparent_edge[ counter_vertices ] = e;
    data_edge_nextparents_of_vertices[ e ][ 1 ] = counter_edges;
    data_edge_nextparents_of_vertices[ counter_edges ][ 0 ] = nullindex;
    
    /* no parent faces yet for the new vertex; will be filled in below */
    
    data_vertex_firstparent_face[ counter_vertices ] = nullindex;
    
    /* no parent faces yet for the front edge; will be filled in below */
    
    data_edge_firstparent_face[ counter_edges ] = nullindex;
    
    
    
    
    
    
    
    
    /* TODO: Bis hierhin geupdated ^^^^^^ */
    
    // run over the list of old faces of the bisected edge 
    // and set the children faces' vertices and parent linkings
    
    for( int of = 0; of < old_faces.size(); of++ ) {
      
      // gather the previous data from the current old face   
      
      int f_old = old_faces[ of ];
      int f_new = counter_faces + of;
      
      int f_e0 = data_face_edges[ f_old ][ 0 ];
      int f_e1 = data_face_edges[ f_old ][ 1 ];
      int f_e2 = data_face_edges[ f_old ][ 2 ];
      
      int f_v0 = data_face_vertices[ f_old ][ 0 ];
      int f_v1 = data_face_vertices[ f_old ][ 1 ];
      int f_v2 = data_face_vertices[ f_old ][ 2 ];
      
      int f_e_n0 = data_face_nextparents_of_edges[ f_old ][ 0 ];
      int f_e_n1 = data_face_nextparents_of_edges[ f_old ][ 1 ];
      int f_e_n2 = data_face_nextparents_of_edges[ f_old ][ 2 ];
      
      int f_v_n0 = data_face_nextparents_of_vertices[ f_old ][ 0 ];
      int f_v_n1 = data_face_nextparents_of_vertices[ f_old ][ 1 ];
      int f_v_n2 = data_face_nextparents_of_vertices[ f_old ][ 2 ];
      
      localindex_of_face_refinementedge[ of ] = ( f_e0 == e ) ? 0 : ( ( f_e1 == e ) ? 1 : 2 );
      assert( data_face_edges[ f_old ][ localindex_of_face_refinementedge[ of ] ] == e );
      
      
      
      if( localindex_of_face_refinementedge[ of ] == 0 ) { // 0 1 
        
        assert( f_v0 == e_back_vertex && f_v1 == e_front_vertex && e == f_e0 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = counter_vertices;
        data_face_vertices[ f_old ][2] = f_v2;
        
        data_face_vertices[ f_new ][0] = counter_vertices;
        data_face_vertices[ f_new ][1] = f_v1;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e;
        data_face_edges[ f_old ][1] = f_e1;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges;
        data_face_edges[ f_new ][1] = counter_edges + 1 + of;
        data_face_edges[ f_new ][2] = f_e2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + of ][1] = f_v2;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][1] = f_v_n1;
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_edges[ f_new ][0] = data_edge_firstparent_face[ counter_edges ];
        data_face_nextparents_of_edges[ f_new ][1] = nullindex;
        data_face_nextparents_of_edges[ f_new ][2] = f_e_n2;
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e2 ] == f_old ) 
          data_edge_firstparent_face[ f_e2 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e2 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v2 ];
        data_vertex_firstparent_edge[ f_v2 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = opposing_vertex_firstparent_edge;
        
        
        
      } else if( localindex_of_face_refinementedge[ of ] == 1 ) { // 0 2 
        
        assert( f_v0 == e_back_vertex && f_v2 == e_front_vertex && e == f_e1 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = f_v1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = f_v1;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = f_e0;
        data_face_edges[ f_old ][1] = e;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = f_e2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = f_v1;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = f_v_n1;
        data_face_nextparents_of_vertices[ f_new ][1] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_edges[ f_new ][0] = nullindex;
        data_face_nextparents_of_edges[ f_new ][1] = f_e_n2;
        data_face_nextparents_of_edges[ f_new ][2] = data_edge_firstparent_face[ counter_edges ];
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e2 ] == f_old ) 
          data_edge_firstparent_face[ f_e2 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e2 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e2 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v1 ];
        data_vertex_firstparent_edge[ f_v1 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = opposing_vertex_firstparent_edge;
        
        
      } else if( localindex_of_face_refinementedge[ of ] == 2 ) { // 1 2 
        
        assert( f_v1 == e_back_vertex && f_v2 == e_front_vertex && e == f_e2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = f_v0;
        data_face_vertices[ f_old ][1] = f_v1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = f_v0;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = f_v2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = f_e0;
        data_face_edges[ f_old ][1] = counter_edges + 1 + of;
        data_face_edges[ f_old ][2] = e;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = f_e1;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = f_v0;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
        
        /* face vertex neighbors */
        data_face_nextparents_of_vertices[ f_old ][0] = counter_faces + of;
        data_face_nextparents_of_vertices[ f_old ][1] = f_v_n1;
        data_face_nextparents_of_vertices[ f_old ][2] = counter_faces + of;
        
        data_face_nextparents_of_vertices[ f_new ][0] = f_v_n0;
        data_face_nextparents_of_vertices[ f_new ][1] = data_vertex_firstparent_face[ counter_vertices ];
        data_face_nextparents_of_vertices[ f_new ][2] = f_v_n2;
        
        data_vertex_firstparent_face[ counter_vertices ] = f_old;
        
        /* face edge neighbors */
        data_face_nextparents_of_edges[ f_old ][0] = f_e_n0;
        data_face_nextparents_of_edges[ f_old ][1] = counter_faces + of;
        data_face_nextparents_of_edges[ f_old ][2] = f_e_n2;
        
        data_face_nextparents_of_edges[ f_new ][0] = nullindex;
        data_face_nextparents_of_edges[ f_new ][1] = f_e_n1;
        data_face_nextparents_of_edges[ f_new ][2] = data_edge_firstparent_face[ counter_edges ];
        
        data_edge_firstparent_face[ counter_edges ] = f_new;
        
        data_edge_firstparent_face[ counter_edges + 1 + of ] = f_old;
        
        /* run over the front outer edge parent face list and replace 'f_old' by 'f_new' */
        if( data_edge_firstparent_face[ f_e1 ] == f_old ) 
          data_edge_firstparent_face[ f_e1 ] = f_new;
        else {
          int current_face = data_edge_firstparent_face[ f_e1 ];
          while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != f_old 
                 &&
                 data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != nullindex )
            current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ];
          assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] != nullindex );
          data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, f_e1 ) ] = f_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ f_v0 ];
        data_vertex_firstparent_edge[ f_v0 ] = counter_edges + 1 + of;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + of ][0] = opposing_vertex_firstparent_edge;
        
      } else {
        
        assert(false);
        
      } 
      
    }
    
    
    
    
    
    
    
    
    
    
    /* Run over the front vertex' parent faces and conduct manipulations */

    int* pointer_to_index = &data_vertex_firstparent_face[ e_front_vertex ];
    
    while( *pointer_to_index != nullindex ) { 
      
      std::vector<int>::iterator it = std::find( old_faces.begin(), old_faces.end(), *pointer_to_index );
      
      if( it != old_faces.end() ) {
        
        assert( *pointer_to_index == *it );
        assert( *it == old_faces[ it - old_faces.begin() ] );
        
        *pointer_to_index = counter_faces + ( it - old_faces.begin() );
        
        std::cout << "manipulate" << std::endl;
        
      } else std::cout << "keep" << std::endl;
      
      int localindex_of_front_vertex = nullindex;
      if( data_face_vertices[ *pointer_to_index ][ 0 ] == e_front_vertex ) localindex_of_front_vertex = 0;
      if( data_face_vertices[ *pointer_to_index ][ 1 ] == e_front_vertex ) localindex_of_front_vertex = 1;
      if( data_face_vertices[ *pointer_to_index ][ 2 ] == e_front_vertex ) localindex_of_front_vertex = 2;
      assert( localindex_of_front_vertex != nullindex );
      
      pointer_to_index = &( data_face_nextparents_of_vertices[ *pointer_to_index ][ localindex_of_front_vertex ] );
      
    }
    
    getcoordinates().append( midcoordinate );
    
    
    
    
    /* TODO: Fill in the data for the case distinctions below */
    
    std::cout << "SECOND BIG LOOP BEGIN" << std::endl;
    
    for( int ot = 0; ot < old_tetrahedra.size(); ot++ ) {
      
      int t_old = old_tetrahedra[ ot ];
      int t_new = counter_tetrahedra + ot;
      
      int f_new = counter_faces + old_faces.size() + ot;
      
      int t_f0 = data_tetrahedron_faces[ t_old ][ 0 ];
      int t_f1 = data_tetrahedron_faces[ t_old ][ 1 ];
      int t_f2 = data_tetrahedron_faces[ t_old ][ 2 ];
      int t_f3 = data_tetrahedron_faces[ t_old ][ 3 ];
      
      int t_e0 = data_tetrahedron_edges[ t_old ][ 0 ];
      int t_e1 = data_tetrahedron_edges[ t_old ][ 1 ];
      int t_e2 = data_tetrahedron_edges[ t_old ][ 2 ];
      int t_e3 = data_tetrahedron_edges[ t_old ][ 3 ];
      int t_e4 = data_tetrahedron_edges[ t_old ][ 4 ];
      int t_e5 = data_tetrahedron_edges[ t_old ][ 5 ];
      
      int t_v0 = data_tetrahedron_vertices[ t_old ][ 0 ];
      int t_v1 = data_tetrahedron_vertices[ t_old ][ 1 ];
      int t_v2 = data_tetrahedron_vertices[ t_old ][ 2 ];
      int t_v3 = data_tetrahedron_vertices[ t_old ][ 3 ];
      
      int t_f_n0 = data_tetrahedron_nextparents_of_faces[ t_old ][ 0 ];
      int t_f_n1 = data_tetrahedron_nextparents_of_faces[ t_old ][ 1 ];
      int t_f_n2 = data_tetrahedron_nextparents_of_faces[ t_old ][ 2 ];
      int t_f_n3 = data_tetrahedron_nextparents_of_faces[ t_old ][ 3 ];
      
      int t_e_n0 = data_tetrahedron_nextparents_of_edges[ t_old ][ 0 ];
      int t_e_n1 = data_tetrahedron_nextparents_of_edges[ t_old ][ 1 ];
      int t_e_n2 = data_tetrahedron_nextparents_of_edges[ t_old ][ 2 ];
      int t_e_n3 = data_tetrahedron_nextparents_of_edges[ t_old ][ 3 ];
      int t_e_n4 = data_tetrahedron_nextparents_of_edges[ t_old ][ 4 ];
      int t_e_n5 = data_tetrahedron_nextparents_of_edges[ t_old ][ 5 ];
      
      int t_v_n0 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 0 ];
      int t_v_n1 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 1 ];
      int t_v_n2 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 2 ];
      int t_v_n3 = data_tetrahedron_nextparents_of_vertices[ t_old ][ 3 ];
      
      localindex_of_tetrahedron_refinementedge[ ot ] = nullindex;
      if( t_e0 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 0;
      if( t_e1 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 1;
      if( t_e2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 2;
      if( t_e3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 3;
      if( t_e4 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 4;
      if( t_e5 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 5;
      assert( localindex_of_tetrahedron_refinementedge[ ot ] != nullindex );
      assert( data_tetrahedron_edges[ t_old ][ localindex_of_tetrahedron_refinementedge[ ot ] ] == e );
      
      
      
      if(        localindex_of_face_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( t_v0 == e_back_vertex && t_v1 == e_front_vertex && e == t_e0 );
        
//         /* Tetrahedron -> Faces */
//         data_tetrahedron_faces[ t_old ][0] = ;
//         data_tetrahedron_faces[ t_old ][1] = ;
//         data_tetrahedron_faces[ t_old ][2] = ;
//         data_tetrahedron_faces[ t_old ][3] = ;
//         
//         data_tetrahedron_faces[ t_new ][0] = ;
//         data_tetrahedron_faces[ t_new ][1] = ;
//         data_tetrahedron_faces[ t_new ][2] = ;
//         data_tetrahedron_faces[ t_new ][3] = ;
//         
//         /* Tetrahedron -> Edges */
//         data_tetrahedron_edges[ t_old ][0] = ;
//         data_tetrahedron_edges[ t_old ][1] = ;
//         data_tetrahedron_edges[ t_old ][2] = ;
//         data_tetrahedron_edges[ t_old ][3] = ;
//         data_tetrahedron_edges[ t_old ][4] = ;
//         data_tetrahedron_edges[ t_old ][5] = ;
//         
//         data_tetrahedron_edges[ t_new ][0] = ;
//         data_tetrahedron_edges[ t_new ][1] = ;
//         data_tetrahedron_edges[ t_new ][2] = ;
//         data_tetrahedron_edges[ t_new ][3] = ;
//         data_tetrahedron_edges[ t_new ][4] = ;
//         data_tetrahedron_edges[ t_new ][5] = ;
//         
//         /* Tetrahedron -> Vertices */
//         data_tetrahedron_vertices[ t_old ][0] = t_v0;
//         data_tetrahedron_vertices[ t_old ][1] = counter_vertices;
//         data_tetrahedron_vertices[ t_old ][2] = t_v2;
//         data_tetrahedron_vertices[ t_old ][3] = t_v3;
//         
//         data_tetrahedron_vertices[ t_new ][0] = counter_vertices;
//         data_tetrahedron_vertices[ t_new ][1] = t_v1;
//         data_tetrahedron_vertices[ t_new ][2] = t_v2;
//         data_tetrahedron_vertices[ t_new ][3] = t_v3;
//         
//         /* Face -> Edges */
//         data_face_edges[ t_old ][0] = ;
//         data_face_edges[ t_old ][1] = ;
//         data_face_edges[ t_old ][2] = ;
//         
//         data_face_edges[ t_new ][0] = ;
//         data_face_edges[ t_new ][1] = ;
//         data_face_edges[ t_new ][2] = ;
//         
//         /* Face -> Vertices */
//         data_face_vertices[ t_old ][0] = ;
//         data_face_vertices[ t_old ][1] = ;
//         data_face_vertices[ t_old ][2] = ;
//         
//         data_face_vertices[ t_new ][0] = ;
//         data_face_vertices[ t_new ][1] = ;
//         data_face_vertices[ t_new ][2] = ;
        
//         /* Neighbors: Tetrahedron -> Faces */
//         /* Neighbors: Tetrahedron -> Edges */
//         /* Neighbors: Tetrahedron -> Vertices */
//         
//         /* Neighbors: Face -> Edges */
//         /* Neighbors: Face -> Vertices */
//         
//         /* run over the front outer face parent tetrahedra and replace 't_old' by 't_new' */
//         if( data_face_firstparent_tetrahedron[ t_f? ] == t_old ) 
//           data_face_firstparent_tetrahedron[ t_f? ] = t_new;
//         else {
//           int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f? ];
//           while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != t_old 
//                  &&
//                  data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != nullindex )
//             current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_e? ) ];
//           assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] != nullindex );
//           data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f? ) ] = t_new;
//         }
//         
//         /* add new parent tetrahedron for the new vertex */
//         int new_vertex_firstparent_tetrahedron = data_vertex_firstparent_tetrahedron[ counter_vertices ];
//         data_vertex_firstparent_tetrahedron[ counter_vertices ] = ;
//         data_edge_nextparents_of_vertices[  ][ ] = new_vertex_firstparent_tetrahedron;
//         
//         /* add new parent face for the opposing edge and its vertices */
        
        
        
      } else if( localindex_of_face_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( t_v0 == e_back_vertex && t_v2 == e_front_vertex && e == t_e1 );
        
      } else if( localindex_of_face_refinementedge[ ot ] == 2 ) { // 0 3 
        
        assert( t_v0 == e_back_vertex && t_v3 == e_front_vertex && e == t_e2 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 3 ) { // 1 2 
        
        assert( t_v1 == e_back_vertex && t_v2 == e_front_vertex && e == t_e3 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 4 ) { // 1 3 
        
        assert( t_v1 == e_back_vertex && t_v3 == e_front_vertex && e == t_e4 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 5 ) { // 2 3 
        
        assert( t_v2 == e_back_vertex && t_v3 == e_front_vertex && e == t_e5 );
    
      } else {
        
        assert( false );
        
      }
       
      
    }
      
    std::cout << "SECOND BIG LOOP END" << std::endl;
    
    
    
    /* SOME FINISHINGS HERE, such as updating the adjacent tetrahedra for front vertex and front edges */
    
        
    
    
    
    
    std::cout << "FINISHED" << std::endl;
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_faces += old_faces.size();
    counter_edges     += 1 + old_faces.size();
    counter_vertices  += 1;
    
    
    
    
    
    /* Done */
    
    check();
    
}

/*
 * * * * * VERTEX BISECTION
 */







