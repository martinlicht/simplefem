
#include <string>
#include <vector>
#include <stack> // TODO: change to something else such as list 
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
    
    
    
    /*********************/
    /*                   */
    /*   GEOMETRY DATA   */
    /*                   */
    /*********************/
    
    getcoordinates().append( midcoordinate );
    
    
    
    
    
    
    /*******************************************/
    /*****                                ******/
    /*****   Assemble sets of simplices   ******/
    /*****   in the micropatch            ******/
    /*****                                ******/
    /*******************************************/
    
    
    /*******************************************/
    /*****                                ******/
    /*****   Assemble sets of simplices   ******/
    /*****   in the MACROpatch            ******/
    /*****                                ******/
    /*******************************************/
    
    
    
    
    
    
    
    
    
    
    /****************************/
    /*                          */
    /*   FILL IN SIMPLEX DATA   */
    /*                          */
    /****************************/
    
    /* vertices of the children of the bisected edge */
    
    data_edge_vertices[ e             ][ 0 ] = e_back_vertex;
    data_edge_vertices[ e             ][ 1 ] = counter_vertices;
      
    data_edge_vertices[ counter_edges ][ 0 ] = counter_vertices;
    data_edge_vertices[ counter_edges ][ 1 ] = e_front_vertex;
    
    
    
    for( int of = 0; of < old_faces.size(); of++ ) {
      
      // gather the previous data from the current old face   
      
      int f_old = old_faces[ of ];
      int f_new = counter_faces + of;
      
      int e_0_1 = data_face_edges[ f_old ][ 0 ];
      int e_0_2 = data_face_edges[ f_old ][ 1 ];
      int e_1_2 = data_face_edges[ f_old ][ 2 ];
      
      int v_0 = data_face_vertices[ f_old ][ 0 ];
      int v_1 = data_face_vertices[ f_old ][ 1 ];
      int v_2 = data_face_vertices[ f_old ][ 2 ];
      
      
      localindex_of_face_refinementedge[ of ] = ( e_0_1 == e ) ? 0 : ( ( e_0_2 == e ) ? 1 : 2 );
      assert( data_face_edges[ f_old ][ localindex_of_face_refinementedge[ of ] ] == e );
      
      
      
      if( localindex_of_face_refinementedge[ of ] == 0 ) { // 0 1 
        
        assert( v_0 == e_back_vertex && v_1 == e_front_vertex && e == e_0_1 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = counter_vertices;
        data_face_vertices[ f_old ][2] = v_2;
        
        data_face_vertices[ f_new ][0] = counter_vertices;
        data_face_vertices[ f_new ][1] = v_1;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e;
        data_face_edges[ f_old ][1] = e_0_2;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges;
        data_face_edges[ f_new ][1] = counter_edges + 1 + of;
        data_face_edges[ f_new ][2] = e_1_2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + of ][1] = v_2;
        
      } else if( localindex_of_face_refinementedge[ of ] == 1 ) { // 0 2 
        
        assert( v_0 == e_back_vertex && v_2 == e_front_vertex && e == e_0_2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = v_1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e_0_1;
        data_face_edges[ f_old ][1] = e;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = e_1_2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = v_1;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
      } else if( localindex_of_face_refinementedge[ of ] == 2 ) { // 1 2 
        
        assert( v_1 == e_back_vertex && v_2 == e_front_vertex && e == e_1_2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = v_1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e_0_1;
        data_face_edges[ f_old ][1] = counter_edges + 1 + of;
        data_face_edges[ f_old ][2] = e;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = e_0_2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = v_0;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
      } else {
        
        assert(false);
        
      } 
      
    }
    
    
    
    for( int ot = 0; ot < old_tetrahedra.size(); ot++ ) {
      
      int t_old = old_tetrahedra[ ot ];
      int t_new = counter_tetrahedra + ot;
      
      int f_new = counter_faces + old_faces.size() + ot;
      
      int f_012 = data_tetrahedron_faces[ t_old ][ 0 ];
      int f_013 = data_tetrahedron_faces[ t_old ][ 1 ];
      int f_023 = data_tetrahedron_faces[ t_old ][ 2 ];
      int f_123 = data_tetrahedron_faces[ t_old ][ 3 ];
      
      int e_0_1 = data_tetrahedron_edges[ t_old ][ 0 ];
      int e_0_2 = data_tetrahedron_edges[ t_old ][ 1 ];
      int e_0_3 = data_tetrahedron_edges[ t_old ][ 2 ];
      int e_1_2 = data_tetrahedron_edges[ t_old ][ 3 ];
      int e_1_3 = data_tetrahedron_edges[ t_old ][ 4 ];
      int e_2_3 = data_tetrahedron_edges[ t_old ][ 5 ];
      
      int v_0 = data_tetrahedron_vertices[ t_old ][ 0 ];
      int v_1 = data_tetrahedron_vertices[ t_old ][ 1 ];
      int v_2 = data_tetrahedron_vertices[ t_old ][ 2 ];
      int v_3 = data_tetrahedron_vertices[ t_old ][ 3 ];
      
      int v_n = counter_vertices; 
      
      localindex_of_tetrahedron_refinementedge[ ot ] = nullindex;
      if( e_0_1 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 0;
      if( e_0_2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 1;
      if( e_0_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 2;
      if( e_1_2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 3;
      if( e_1_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 4;
      if( e_2_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 5;
      assert( localindex_of_tetrahedron_refinementedge[ ot ] != nullindex );
      assert( data_tetrahedron_edges[ t_old ][ localindex_of_tetrahedron_refinementedge[ ot ] ] == e );
      
      
      
      if(        localindex_of_face_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( v_0 == e_back_vertex && v_1 == e_front_vertex && e == e_0_1 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        
        int e_0_n = e_0_1;
        int e_n_1 = v_n;
        int e_n_2 = counter_edges + 1 + index_f_012;
        int e_n_3 = counter_edges + 1 + index_f_013;
        
        int f_0n2 = f_012;
        int f_n12 = counter_faces + index_f_012;
        int f_0n3 = f_013;
        int f_n13 = counter_faces + index_f_013;
        int f_n23 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_0n2;
        data_tetrahedron_faces[ t_old ][1] = f_0n3;
        data_tetrahedron_faces[ t_old ][2] = f_023;
        data_tetrahedron_faces[ t_old ][3] = f_n23;
        
        data_tetrahedron_faces[ t_new ][0] = f_n12;
        data_tetrahedron_faces[ t_new ][1] = f_n13;
        data_tetrahedron_faces[ t_new ][2] = f_n23;
        data_tetrahedron_faces[ t_new ][3] = f_123;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_n;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_n_2;
        data_tetrahedron_edges[ t_old ][4] = e_n_3;
        data_tetrahedron_edges[ t_old ][5] = e_2_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_n_1;
        data_tetrahedron_edges[ t_new ][1] = e_n_2;
        data_tetrahedron_edges[ t_new ][2] = e_n_3;
        data_tetrahedron_edges[ t_new ][3] = e_1_2;
        data_tetrahedron_edges[ t_new ][4] = e_1_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_n;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_n;
        data_tetrahedron_vertices[ t_new ][1] = v_1;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_n_2;
        data_face_edges[ f_new ][1] = e_n_3;
        data_face_edges[ f_new ][2] = e_2_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_n;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_3;
        
      } else if( localindex_of_face_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( v_0 == e_back_vertex && v_2 == e_front_vertex && e == e_0_2 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        
        int e_0_n = e_0_2;
        int e_n_2 = v_n;
        int e_1_n = counter_edges + 1 + index_f_012;
        int e_n_3 = counter_edges + 1 + index_f_023;
        
        int f_01n = f_012;
        int f_1n2 = counter_faces + index_f_012;
        int f_0n3 = f_023;
        int f_2n3 = counter_faces + index_f_023;
        int f_1n3 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_01n;
        data_tetrahedron_faces[ t_old ][1] = f_013;
        data_tetrahedron_faces[ t_old ][2] = f_0n3;
        data_tetrahedron_faces[ t_old ][3] = f_1n3;
        
        data_tetrahedron_faces[ t_new ][0] = f_1n2;
        data_tetrahedron_faces[ t_new ][1] = f_1n3;
        data_tetrahedron_faces[ t_new ][2] = f_123;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_n;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_1_n;
        data_tetrahedron_edges[ t_old ][4] = e_1_3;
        data_tetrahedron_edges[ t_old ][5] = e_n_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_1_n;
        data_tetrahedron_edges[ t_new ][1] = e_1_2;
        data_tetrahedron_edges[ t_new ][2] = e_1_3;
        data_tetrahedron_edges[ t_new ][3] = e_n_2;
        data_tetrahedron_edges[ t_new ][4] = e_n_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_n;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_1;
        data_tetrahedron_vertices[ t_new ][1] = v_n;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_1_n;
        data_face_edges[ f_new ][1] = e_1_3;
        data_face_edges[ f_new ][2] = e_n_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = v_n;
        data_face_vertices[ f_new ][2] = v_3;
        
      } else if( localindex_of_face_refinementedge[ ot ] == 2 ) { // 0 3 
        
        assert( v_0 == e_back_vertex && v_3 == e_front_vertex && e == e_0_3 );
        
        /* prepare stuff */ 
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        
        int e_0_n = e_0_3;
        int e_n_3 = v_n;
        int e_1_n = counter_edges + 1 + index_f_013;
        int e_2_n = counter_edges + 1 + index_f_023;
        
        int f_01n = f_013;
        int f_1n3 = counter_faces + index_f_013;
        int f_02n = f_023;
        int f_2n3 = counter_faces + index_f_023;
        int f_12n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_12n;
        data_tetrahedron_faces[ t_new ][1] = f_123;
        data_tetrahedron_faces[ t_new ][2] = f_1n3;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_1_2;
        data_tetrahedron_edges[ t_new ][1] = e_1_n;
        data_tetrahedron_edges[ t_new ][2] = e_1_3;
        data_tetrahedron_edges[ t_new ][3] = e_2_n;
        data_tetrahedron_edges[ t_new ][4] = e_2_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_1;
        data_tetrahedron_vertices[ t_new ][1] = v_2;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_1_2;
        data_face_edges[ f_new ][1] = e_1_n;
        data_face_edges[ f_new ][2] = e_2_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else if( localindex_of_face_refinementedge[ ot ] == 3 ) { // 1 2 
        
        assert( v_1 == e_back_vertex && v_2 == e_front_vertex && e == e_1_2 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        int e_1_n = e_1_2;
        int e_n_2 = v_n;
        int e_0_n = counter_edges + 0 + index_f_012;
        int e_n_3 = counter_edges + 0 + index_f_123;
        
        int f_01n = f_012;
        int f_0n2 = counter_faces + index_f_012;
        int f_1n3 = f_123;
        int f_n23 = counter_faces + index_f_123;
        int f_0n3 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_01n;
        data_tetrahedron_faces[ t_old ][1] = f_013;
        data_tetrahedron_faces[ t_old ][2] = f_0n3;
        data_tetrahedron_faces[ t_old ][3] = f_1n3;
        
        data_tetrahedron_faces[ t_new ][0] = f_0n2;
        data_tetrahedron_faces[ t_new ][1] = f_0n3;
        data_tetrahedron_faces[ t_new ][2] = f_023;
        data_tetrahedron_faces[ t_new ][3] = f_n23;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_n;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_1_n;
        data_tetrahedron_edges[ t_old ][4] = e_1_3;
        data_tetrahedron_edges[ t_old ][5] = e_n_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_n;
        data_tetrahedron_edges[ t_new ][1] = e_0_2;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_n_2;
        data_tetrahedron_edges[ t_new ][4] = e_n_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_n;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_n;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_n;
        data_face_edges[ f_new ][1] = e_0_3;
        data_face_edges[ f_new ][2] = e_n_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_n;
        data_face_vertices[ f_new ][2] = v_3;
    
      } else if( localindex_of_face_refinementedge[ ot ] == 4 ) { // 1 3 
        
        assert( v_1 == e_back_vertex && v_3 == e_front_vertex && e == e_1_3 );
        
        /* prepare stuff */ 
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        int e_1_n = e_1_3;
        int e_n_3 = v_n;
        int e_0_n = counter_edges + 0 + index_f_013;
        int e_2_n = counter_edges + 0 + index_f_123;
        
        int f_01n = f_013;
        int f_0n3 = counter_faces + index_f_013;
        int f_12n = f_123;
        int f_2n3 = counter_faces + index_f_123;
        int f_02n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_02n;
        data_tetrahedron_faces[ t_new ][1] = f_023;
        data_tetrahedron_faces[ t_new ][2] = f_0n3;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_2;
        data_tetrahedron_edges[ t_new ][1] = e_0_n;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_2_n;
        data_tetrahedron_edges[ t_new ][4] = e_2_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_2;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_2;
        data_face_edges[ f_new ][1] = e_0_n;
        data_face_edges[ f_new ][2] = e_2_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else if( localindex_of_face_refinementedge[ ot ] == 5 ) { // 2 3 
        
        assert( v_2 == e_back_vertex && v_3 == e_front_vertex && e == e_2_3 );
        
        /* prepare stuff */ 
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        int e_2_n = e_2_3;
        int e_n_3 = v_n;
        int e_0_n = counter_edges + 0 + index_f_023;
        int e_1_n = counter_edges + 0 + index_f_123;
        
        int f_02n = f_023;
        int f_0n3 = counter_faces + index_f_023;
        int f_12n = f_123;
        int f_1n3 = counter_faces + index_f_123;
        int f_01n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_01n;
        data_tetrahedron_faces[ t_new ][1] = f_013;
        data_tetrahedron_faces[ t_new ][2] = f_0n3;
        data_tetrahedron_faces[ t_new ][3] = f_1n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_1;
        data_tetrahedron_edges[ t_new ][1] = e_0_n;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_1_n;
        data_tetrahedron_edges[ t_new ][4] = e_1_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_1;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_1;
        data_face_edges[ f_new ][1] = e_0_n;
        data_face_edges[ f_new ][2] = e_1_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_1;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else {
        
        assert( false );
        
      }
       
      
    }
      
    
    
    
    
    /***********************************************/
    /*****                                    ******/
    /*****   Delete the adjaceny links        ******/
    /*****   in the macropatch                ******/
    /*****   regarding micropatch simplices   ******/
    /*****                                    ******/
    /***********************************************/
    
    
    /***********************************************/
    /*****                                    ******/
    /*****   Rebuild the adjaceny links       ******/
    /*****   in the macropatch                ******/
    /*****   regarding micropatch simplices   ******/
    /*****                                    ******/
    /***********************************************/
    
    
    
    
    
    /**************************/
    /*     ADJACENCY LINKS    */
    /**************************/
    
    // Run over the simplices in the macropatch 
    // and reassemble the links corresponding 
    // to simplices in the micropatch
    // or newly created. 
    
    
    /**************************************************/
    
    
    
        
    
    
    
    
    std::cout << "FINISHED" << std::endl;
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_tetrahedra += old_tetrahedra.size();
    counter_faces      += old_faces.size() + old_tetrahedra.size();
    counter_edges      += 1 + old_faces.size();
    counter_vertices   += 1;
    
    
    
    
    
    /* Done */
    
    check();
    
}

/*
 * * * * * VERTEX BISECTION
 */





void MeshSimplicial3D::longest_edge_bisection( std::vector<int> edges )
{
    check();
    
    const int old_vertex_count = counter_vertices;
    
    
    /* 0. check the input */
    
    for( int& e : edges )
        assert( 0 <= edges[e] && edges[e] < counter_edges );
    
    
    /* 1. create stack for the edges to be bisected, and fill in first batch */
    
    std::stack<int> todostack;
    
    for( int& e : edges )
        todostack.push( e ); // put e on top either by inserting or pulling it up!
        
    
    /* 2. conduct the main loop of the refinement algorithm */
    
    while( ! todostack.empty() )
    {
        
        // as long the stack is not empty,
        // pick the top edge and make the following case distinction
        // a) the edge index belongs to an edge already bisected and can be ignored 
        // b) if the top edge is longer than its neighbors, bisect and pop
        // c) else, push the longest edge of each parent simplex
        
        int e = todostack.top();
        
        // to check whether e belongs to an edge that has already been bisected,
        // we check whether one of the vertices belongs to the new vertices 
        
        // TODO: This is actually unsafe
        // A better solution is to replace the stack by a linked list 
        // and emulate the stack behavior yourself. 
        // Then multiple occurences of the same edge 
        // can be removed by an STL algorithm 
        
        
        if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count ) {
            
            todostack.pop();
        
        } else {
            
            Float length_e = get_edge_length( e );
                
            // run over neighbor tetrahedra and check for longer edges 
            
            for(
                int t = get_edge_firstparent_tetrahedron( e );
                t != nullindex; 
                t = get_edge_nextparent_tetrahedron( e, t )
            )
            for( int ei = 0; ei < 6; ei++ )
                if( e != get_tetrahedron_edge( t, ei ) && get_edge_length( get_tetrahedron_edge(t,ei) ) > length_e )
                    todostack.push( get_tetrahedron_edge( t, ei ) );
            
            // if top edge is still the same, bisect 
            
            if( e == todostack.top() )
            {
                todostack.pop();
                bisect_edge( e );
            }
            
        }
        
    }
    
    
    // fiinished!
    
    check();
}









