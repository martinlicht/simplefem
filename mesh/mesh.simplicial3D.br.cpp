
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
        
      } else if( localindex_of_face_refinementedge[ ot ] == 2 ) { // 0 3 
        
        assert( v_0 == e_back_vertex && v_3 == e_front_vertex && e == e_0_3 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 3 ) { // 1 2 
        
        assert( v_1 == e_back_vertex && v_2 == e_front_vertex && e == e_1_2 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 4 ) { // 1 3 
        
        assert( v_1 == e_back_vertex && v_3 == e_front_vertex && e == e_1_3 );
    
      } else if( localindex_of_face_refinementedge[ ot ] == 5 ) { // 2 3 
        
        assert( v_2 == e_back_vertex && v_3 == e_front_vertex && e == e_2_3 );
    
      } else {
        
        assert( false );
        
      }
       
      
    }
      
    
    
    
    
    
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
    
    counter_faces += old_faces.size();
    counter_edges     += 1 + old_faces.size();
    counter_vertices  += 1;
    
    
    
    
    
    /* Done */
    
    check();
    
}

/*
 * * * * * VERTEX BISECTION
 */








