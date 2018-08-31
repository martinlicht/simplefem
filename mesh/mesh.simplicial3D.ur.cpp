
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
 *
 *  Split Tetrahedra: 8T 
 *  
 *    00 11 22 33
 *       00 01 02 03
 *       01 11 12 13
 *       02 12 22 23
 *       03 13 23 33
 *       
 *       01 02 03 13
 *       01 02 12 13 
 *       02 03 13 23
 *       02 12 13 23
 *
 *  Split Faces:    4F
 *    
 *    00 11 22
 *       00 01 02
 *       01 11 12
 *       02 12 22
 *       01 02 12
 *    00 11 33
 *       00 01 03
 *       01 11 13
 *       03 13 33
 *       01 03 13
 *    00 22 33
 *       00 02 03
 *       02 22 23
 *       03 23 33
 *       02 03 23
 *    11 22 33
 *       11 12 13
 *       12 22 23
 *       13 23 33
 *       12 13 23
 *       
 *  Completely New Faces: 8T
 *       
 *       01 02 03
 *       01 12 13
 *       02 12 23
 *       03 13 23
 *       
 *       01 02 13
 *       02 03 13
 *       02 12 13
 *       02 13 23
 *       
 * 
 *  Split Edges: 2 E 
 *    
 *    00 11
 *       00 01
 *       01 11
 *    00 22
 *       00 02
 *       02 22
 *    00 33
 *       00 03
 *       03 33
 *    11 22
 *       11 12
 *       12 22
 *    11 33
 *       11 13
 *       13 33
 *    22 33
 *       22 23
 *       23 33
 *    
 *  Completely New Edge: 3F + 1T
 *       
 *       
 *       01 02
 *       01 12
 *       02 12
 *       
 *       01 03
 *       01 13
 *       03 13
 *       
 *       02 03
 *       02 23
 *       03 23
 *       
 *       12 13
 *       12 23
 *       13 23
 *       
 * 
 * 
 *       02 13
 *      
 *
 *
 *
 */




// TODO: Debug the uniform refinement method 

void MeshSimplicial3D::uniformrefinement()
{
    check();
    
    
    int new_counter_tetrahedra = 8 * counter_tetrahedra;
    int new_counter_faces      = 4 * counter_faces + 8 * counter_tetrahedra;
    int new_counter_edges      = 2 * counter_edges + 3  * counter_faces + 1 * counter_tetrahedra;
    int new_counter_vertices   = 1 * counter_vertices + 1 * counter_edges;
    
    
    /* resize the arrays */
    
    /* tetrahedron -> face */
    
    data_tetrahedron_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_firstparent_tetrahedron.resize( new_counter_faces, nullindex );
    
    data_tetrahedron_nextparents_of_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    
    /* tetrahedron -> edge */
    
    data_tetrahedron_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_tetrahedron.resize( new_counter_edges, nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    
    /* tetrahedron -> vertex */
    
    data_tetrahedron_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_tetrahedron.resize( new_counter_vertices, nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    
    /* face -> edge */
    
    data_face_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_face.resize( new_counter_edges, nullindex );
    
    data_face_nextparents_of_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* face -> vertex */
    
    data_face_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_face.resize( new_counter_vertices, nullindex );
    
    data_face_nextparents_of_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* edge -> vertex */
    
    data_edge_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    data_vertex_firstparent_edge.resize( new_counter_vertices, nullindex );
    
    data_edge_nextparents_of_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    
    /* coordinates */
    
    getcoordinates().addcoordinates( counter_edges );
    
    
    
    
    
    
    /* 0. create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    
    
    
    /********************************/
    /***   VERTICES AND EDGES    ****/
    /********************************/
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO EDGES ****/
    
    /* 1. for each old vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex && 0 <= p && p < counter_edges );
      
      int vi = data_edge_vertices[p][0] == v ? 0 : 1;
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent_edge[v] = vi * counter_edges + p;
    }
    
    
    /* 2. for each old edge, relocate the data of the old vertices' old parent edges */
    
    for( int e  = 0; e  < counter_edges;  e++ )
    for( int vi = 0; vi <=            1; vi++ )
    {
      int q = data_edge_nextparents_of_vertices[e][vi];
      
      int v = data_edge_vertices[e][vi];
      
      assert( v != nullindex && 0 <= v && v < counter_vertices );
      
      if( q == nullindex ) {
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = nullindex;
        
      } else {
        
        assert( 0 <= q && q < counter_edges );
        
        int vinp = data_edge_vertices[q][0] == v ? 0 : 1;
        
        assert( data_edge_vertices[q][0] == v || data_edge_vertices[q][1] == v );
        assert( data_edge_vertices[q][vinp] == v );
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = vinp * counter_edges + q;
      
      } 
      
    }
    
    
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO EDGES ****/
    
    /* 1.
     * for each new vertex (which is in the middle of an old edge), 
     * set the first and second parent edge from the old edge 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent_edge[counter_vertices + e] = e;
      
      data_edge_nextparents_of_vertices[ 0 * counter_edges + e ][1] = counter_edges + e;
      data_edge_nextparents_of_vertices[ 1 * counter_edges + e ][0] = nullindex;
    }
    
    
    
    /* 2.
     * for each old edge, run over the adjacent faces 
     * and add the corresponding new edges to the list of 
     * parent edges of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int f = data_edge_firstparent_face[e];
      
      while( f != nullindex ) {
        
        int ei   = nullindex;
        int e_1  = nullindex; 
        int e_2  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        
        if(        data_face_edges[f][0] == e ) {
          ei = 0; e_1 = 0; e_2 = 1; vi_1 = 0; vi_2 = 0; 
        } else if( data_face_edges[f][1] == e ) {
          ei = 1; e_1 = 0; e_2 = 2; vi_1 = 1; vi_2 = 0; 
        } else if( data_face_edges[f][2] == e ) {
          ei = 2; e_1 = 1; e_2 = 2; vi_1 = 1; vi_2 = 1; 
        } else
          assert(false);
        
        assert( ei  != nullindex && e_1 != nullindex && e_2 != nullindex );
        
        int old_first_parent = data_vertex_firstparent_edge[ counter_vertices + e ];
        
        data_vertex_firstparent_edge[ counter_vertices + e ]
          = 2 * counter_edges + e_1 * counter_faces + f;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_1 * counter_faces + f ][ vi_1 ]
          = 2 * counter_edges + e_2 * counter_faces + f;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_2 * counter_faces + f ][ vi_2 ]
          = old_first_parent;
        
        f = data_face_nextparents_of_edges[ f ][ ei ];
        
      }
      
    }
    
    
    
    
    /* 3.
     * for each tetrahedron, include the single interior edge 02 13
     * in the parent lists of its two vertices
     */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      /* 01 02 03 12 13 23 */
      
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      
      int ne = 2 * counter_edges + 3 * counter_faces + t;
      
      int fp_v02 = data_vertex_firstparent_edge[ v02 ];
      int fp_v13 = data_vertex_firstparent_edge[ v13 ];
      
      data_vertex_firstparent_edge[ v02 ] = ne;
      data_vertex_firstparent_edge[ v13 ] = ne;
      
      data_edge_nextparents_of_vertices[ ne ][ 0 ] = fp_v02;
      data_edge_nextparents_of_vertices[ ne ][ 1 ] = fp_v13;
      
    }
    
    
    
    
    /*** SET THE VERTICES OF ALL EDGES ****/
    
    /* 1. for each edge created from an old edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    /* 2. for each face, set the vertices of the new edges inside the face */
    
    for( int f = 0; f < counter_faces; f++ )
    {
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][1];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][1];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
    }
    
    /* 3. for each tetrahedron, set the internal vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      /* 01 02 03 12 13 23 */
      
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int ne = 2 * counter_edges + 3 * counter_faces + t;
      
      data_edge_vertices[ ne ][ 0 ] = v02;
      data_edge_vertices[ ne ][ 1 ] = v13;
      
    }
    
    
    
    
    
    /****************************************/
    /***   VERTICES AND EDGES AND FACES  ****/
    /****************************************/
    
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO FACES ****/
    
    /* 1. for each old vertex, set the new parent face */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_face[v];
      
      assert( p != nullindex && 0 <= p && p < counter_faces );
      
      int vi = data_face_vertices[p][0] == v ? 0 : data_face_vertices[p][1] == v ? 1 : 2;
      
      assert( data_face_vertices[p][0] == v || data_face_vertices[p][1] == v || data_face_vertices[p][2] == v );
      assert( data_face_vertices[p][vi] == v );
      
      data_vertex_firstparent_face[v] = vi * counter_faces + p;
    }
    
    
    /* 2. for each old face, relocate the data of the old vertices' parent face */
    
    for( int f  = 0; f  < counter_faces;  f++ )
    for( int vi = 0; vi <             3; vi++ )
    {
      
      int q = data_face_nextparents_of_vertices[f][vi];
      
      int v = data_face_vertices[f][vi];
      
      if( q == nullindex ) {
        
        data_face_nextparents_of_vertices[ vi * counter_faces + f ][vi] = nullindex;
        
      } else {
        
        int vinp = data_face_vertices[q][0] == v ? 0 : data_face_vertices[q][1] == v ? 1 : 2;
        
        assert( data_face_vertices[q][0] == v || data_face_vertices[q][1] == v || data_face_vertices[q][2] == v );
        assert( data_face_vertices[q][vinp] == v );
        
        data_face_nextparents_of_vertices[ vi * counter_faces + f ][vi] = vinp * counter_faces + q;
      
      } 
      
    }
    
    
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO FACES ****/
    
    /* 1.
     * for each old edge, run over the adjacent old faces 
     * and add the corresponding new faces to the list of 
     * parent faces of new vertex.
     */
   
   /* TODO: the following code is currently considered legacy */
    
//     for( int e = 0; e < counter_edges; e++ )
//     {
//       int f = data_edge_firstparent_face[e];
//       
//       while( f != nullindex ) {
//         
//         int ei   = nullindex;
//         int f_1  = nullindex; 
//         int f_2  = nullindex;
//         int f_3  = nullindex;
//         int vi_1 = nullindex;
//         int vi_2 = nullindex;
//         int vi_3 = nullindex;
//         
//         // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
//         // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
//         
//         if(        data_face_edges[f][0] == e ) {
//           ei = 0; 
//           f_1 = 0; f_2 = 3; f_3 = 1; vi_1 = 1; vi_2 = 0; vi_3 = 0;  
//         } else if( data_face_edges[f][1] == e ) {
//           ei = 1; 
//           f_1 = 0; f_2 = 3; f_3 = 2; vi_1 = 2; vi_2 = 1; vi_3 = 0;  
//         } else if( data_face_edges[f][2] == e ) {
//           ei = 2; 
//           f_1 = 1; f_2 = 3; f_3 = 2; vi_1 = 2; vi_2 = 2; vi_3 = 1; 
//         } else
//           assert(false);
//         
//         int old_first_parent = data_vertex_firstparent_face[ counter_vertices + e ];
//         
//         data_vertex_firstparent_face[ counter_vertices + e ]
//           = f_1 * counter_faces + f;
//         
//         if( f_1 != 0 ) assert( data_face_nextparents_of_vertices[ f_1 * counter_faces + f ][ vi_1 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_1 * counter_faces + f ][ vi_1 ]
//           = f_2 * counter_faces + f;
//         
//         if( f_2 != 0 ) assert( data_face_nextparents_of_vertices[ f_2 * counter_faces + f ][ vi_2 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_2 * counter_faces + f ][ vi_2 ]
//           = f_3 * counter_faces + f;
//         
//         if( f_3 != 0 ) assert( data_face_nextparents_of_vertices[ f_3 * counter_faces + f ][ vi_3 ] == nullindex );
//         data_face_nextparents_of_vertices[ f_3 * counter_faces + f ][ vi_3 ]
//           = old_first_parent;
//         
//         f = data_face_nextparents_of_edges[ f ][ ei ];
//         
//       }
//       
//     }
    
    for( int f = 0; f < counter_faces; f++ )
    {
        
        int e_01 = data_face_edges[ f ][0];
        int e_02 = data_face_edges[ f ][1];
        int e_12 = data_face_edges[ f ][2];
        
        int v01 = counter_vertices + e_01;
        int v02 = counter_vertices + e_02;
        int v12 = counter_vertices + e_12;
        
        int f_00_01_02 = 0 * counter_faces + f;
        int f_01_11_12 = 1 * counter_faces + f;
        int f_02_12_22 = 2 * counter_faces + f;
        int f_01_02_12 = 3 * counter_faces + f;
        
        int fp_v01 = data_vertex_firstparent_face[ v01 ];
        int fp_v02 = data_vertex_firstparent_face[ v02 ];
        int fp_v12 = data_vertex_firstparent_face[ v12 ];
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
        data_vertex_firstparent_face[ v01 ] = f_00_01_02;
        data_face_nextparents_of_vertices[ f_00_01_02 ][ 1 ] = f_01_11_12;
        data_face_nextparents_of_vertices[ f_01_11_12 ][ 0 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 0 ] = fp_v01;
        
        data_vertex_firstparent_face[ v02 ] = f_00_01_02;
        data_face_nextparents_of_vertices[ f_00_01_02 ][ 2 ] = f_02_12_22;
        data_face_nextparents_of_vertices[ f_02_12_22 ][ 0 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 1 ] = fp_v02;
        
        data_vertex_firstparent_face[ v12 ] = f_01_11_12;
        data_face_nextparents_of_vertices[ f_01_11_12 ][ 2 ] = f_02_12_22;
        data_face_nextparents_of_vertices[ f_02_12_22 ][ 1 ] = f_01_02_12;
        data_face_nextparents_of_vertices[ f_01_02_12 ][ 2 ] = fp_v12;
        
    }
    
    
    
    
    
    
    /* for each tetrahedron, add the new interior faces to the parent lists of the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
       /*
        *       01 02 03
        *       01 12 13
        *       02 12 23
        *       03 13 23
        *       
        *       01 02 13
        *       02 03 13
        *       02 12 13
        *       02 13 23
        */
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      int fp_v01 = data_vertex_firstparent_face[ v01 ];
      int fp_v02 = data_vertex_firstparent_face[ v02 ];
      int fp_v03 = data_vertex_firstparent_face[ v03 ];
      int fp_v12 = data_vertex_firstparent_face[ v12 ];
      int fp_v13 = data_vertex_firstparent_face[ v13 ];
      int fp_v23 = data_vertex_firstparent_face[ v23 ];
      
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      data_vertex_firstparent_face[ v01 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 0 ] = f_01_12_13; 
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_vertices[ f_01_02_13 ][ 0 ] = fp_v01; 
      
      data_vertex_firstparent_face[ v02 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 1 ] = f_02_12_23; 
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_vertices[ f_01_02_13 ][ 1 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 0 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 0 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 0 ] = fp_v02; 
      
      data_vertex_firstparent_face[ v03 ] = f_01_02_03;
      data_face_nextparents_of_vertices[ f_01_02_03 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 0 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 1 ] = fp_v03; 
      
      data_vertex_firstparent_face[ v12 ] = f_01_12_13;
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 1 ] = f_02_12_23; 
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 1 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 1 ] = fp_v12; 
      
      data_vertex_firstparent_face[ v13 ] = f_01_12_13;
      data_face_nextparents_of_vertices[ f_01_12_13 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 1 ] = f_02_03_13; 
      data_face_nextparents_of_vertices[ f_02_03_13 ][ 2 ] = f_02_12_13; 
      data_face_nextparents_of_vertices[ f_02_12_13 ][ 2 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 1 ] = fp_v13; 
      
      data_vertex_firstparent_face[ v23 ] = f_02_12_23;
      data_face_nextparents_of_vertices[ f_02_12_23 ][ 2 ] = f_03_13_23; 
      data_face_nextparents_of_vertices[ f_03_13_23 ][ 2 ] = f_02_13_23; 
      data_face_nextparents_of_vertices[ f_02_13_23 ][ 2 ] = fp_v23; 
      
    }
    
    
    
    /*** TREAT THE BISECTED EDGES AND THEIR CONNECTION TO FACES ****/
    
    
    /* for each bisected edge, set the new first parent faces of the children edges */
    for( int e = 0; e < counter_edges; e++ )
    {
      int p = data_edge_firstparent_face[e];
      
      assert( p != nullindex );
      
      int ei        = nullindex;
      int nfp_back  = nullindex;
      int nfp_front = nullindex;
      
      if( data_face_edges[p][0] == e ){
        ei = 0; nfp_back = 0; nfp_front = 1;
      } else if( data_face_edges[p][1] == e ) {
        ei = 1; nfp_back = 0; nfp_front = 2;
      } else if( data_face_edges[p][2] == e ) {
        ei = 2; nfp_back = 1; nfp_front = 2;
      } else 
        assert(false);
      
      assert( ei != nullindex );
      assert( data_face_edges[p][0] == e || data_face_edges[p][1] == e || data_face_edges[p][2] == e );
      assert( data_face_edges[p][ei] == e );
      
      data_edge_firstparent_face[ 0 * counter_edges + e ] = nfp_back  * counter_faces + p;
      data_edge_firstparent_face[ 1 * counter_edges + e ] = nfp_front * counter_faces + p;
      
    }
    
    
    /* for each face, relocate the data of the old edges' parent face */
    /* additionally, set the new parents */
    for( int f  = 0; f  < counter_faces;  f++ )
    for( int ei = 0; ei <             3; ei++ )
    {
      int e = data_face_edges[f][ei];
      
      int f_back  = nullindex;
      int f_front = nullindex;
      int e_back  = nullindex;
      int e_front = nullindex;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      
      if( ei == 0 ){
        f_back = 0; f_front = 1; e_back = 0; e_front = 0; 
      } else if( ei == 1 ) {
        f_back = 0; f_front = 2; e_back = 1; e_front = 1; 
      } else if( ei == 2 ) {
        f_back = 1; f_front = 2; e_back = 2; e_front = 2; 
      } else 
        assert(false);
      
      
      int q = data_face_nextparents_of_edges[f][ei];
      
      if( q == nullindex ) {
        
        data_face_nextparents_of_edges[ f_back  * counter_faces + f ][ e_back  ] = nullindex;
        data_face_nextparents_of_edges[ f_front * counter_faces + f ][ e_front ] = nullindex;
        
      } else if( q != nullindex ) {
        
        int q_ei        = nullindex;
        int q_nfp_back  = nullindex;
        int q_nfp_front = nullindex;
        
        if(        data_face_edges[q][0] == e ){
          q_ei = 0; q_nfp_back = 0; q_nfp_front = 1;
        } else if( data_face_edges[q][1] == e ) {
          q_ei = 1; q_nfp_back = 0; q_nfp_front = 2;
        } else if( data_face_edges[q][2] == e ) {
          q_ei = 2; q_nfp_back = 1; q_nfp_front = 2;
        } else 
          assert(false);
        
        assert( q_ei != nullindex );
        assert( data_face_edges[q][0] == e || data_face_edges[q][1] == e || data_face_edges[q][2] == e );
        assert( data_face_edges[q][q_ei] == e );
        
        data_face_nextparents_of_edges[ f_back  * counter_faces + f ][ e_back  ] = q_nfp_back  * counter_faces + q;
        data_face_nextparents_of_edges[ f_front * counter_faces + f ][ e_front ] = q_nfp_front * counter_faces + q;
        
      } 
      
    }
    
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO SPLIT FACES ****/
    
    /* for each face, run over the new edges and add firstparents and parents */
    for( int f = 0; f < counter_faces; f++ )
    {
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      // new edges [ 01 02 ], [ 01 12 ], [ 02 12 ]
      
      data_edge_firstparent_face[ 2 * counter_edges + 0 * counter_faces + f ] = 3 * counter_faces + f;
      data_edge_firstparent_face[ 2 * counter_edges + 1 * counter_faces + f ] = 3 * counter_faces + f;
      data_edge_firstparent_face[ 2 * counter_edges + 2 * counter_faces + f ] = 3 * counter_faces + f;
      
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][0] = 0 * counter_faces + f;
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][1] = 1 * counter_faces + f;
      data_face_nextparents_of_edges[ 3 * counter_faces + f ][2] = 2 * counter_faces + f;
      
      data_face_nextparents_of_edges[ 0 * counter_faces + f ][2] = nullindex;
      data_face_nextparents_of_edges[ 1 * counter_faces + f ][1] = nullindex;
      data_face_nextparents_of_edges[ 2 * counter_faces + f ][0] = nullindex;
      
    }
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO INTERIOR FACES ****/
    /*** TREAT THE SINGLE INTERIOR EDGE AND ITS CONNECTION TO INTERIOR FACES ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
       /*
        *       01 02 03
        *       01 12 13
        *       02 12 23
        *       03 13 23
        *       
        *       01 02 13
        *       02 03 13
        *       02 12 13
        *       02 13 23
        */
      
      int f_012 = data_tetrahedron_faces[t][0];
      int f_013 = data_tetrahedron_faces[t][1];
      int f_023 = data_tetrahedron_faces[t][2];
      int f_123 = data_tetrahedron_faces[t][3];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
      
      int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
      int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
      int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
      
      int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
      int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
      int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
      
      int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
      int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
      int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
      
      int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
      
      
      int fp_e_01_02 = data_edge_firstparent_face[ e_01_02 ];
      int fp_e_01_12 = data_edge_firstparent_face[ e_01_12 ];
      int fp_e_02_12 = data_edge_firstparent_face[ e_02_12 ];
      
      int fp_e_01_03 = data_edge_firstparent_face[ e_01_03 ];
      int fp_e_01_13 = data_edge_firstparent_face[ e_01_13 ];
      int fp_e_03_13 = data_edge_firstparent_face[ e_03_13 ];
      
      int fp_e_02_03 = data_edge_firstparent_face[ e_02_03 ];
      int fp_e_02_23 = data_edge_firstparent_face[ e_02_23 ];
      int fp_e_03_23 = data_edge_firstparent_face[ e_03_23 ];
      
      int fp_e_12_13 = data_edge_firstparent_face[ e_12_13 ];
      int fp_e_12_23 = data_edge_firstparent_face[ e_12_23 ];
      int fp_e_13_23 = data_edge_firstparent_face[ e_13_23 ];
      
      int fp_e_02_13 = data_edge_firstparent_face[ e_02_13 ];
      
      
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      data_edge_firstparent_face[ e_01_02 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 0 ] = f_01_02_13; 
      data_face_nextparents_of_edges[ f_01_02_13 ][ 0 ] = fp_e_01_02; 
      
      data_edge_firstparent_face[ e_01_12 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 0 ] = fp_e_01_12;
      
      data_edge_firstparent_face[ e_02_12 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 0 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 0 ] = fp_e_02_12;
      
      
      data_edge_firstparent_face[ e_01_03 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 1 ] = fp_e_01_03;
      
      data_edge_firstparent_face[ e_01_13 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 1 ] = f_01_02_13;
      data_face_nextparents_of_edges[ f_01_02_13 ][ 1 ] = fp_e_01_13;
      
      data_edge_firstparent_face[ e_03_13 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 0 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 2 ] = fp_e_03_13;
      
      
      data_edge_firstparent_face[ e_02_03 ] = f_01_02_03;
      data_face_nextparents_of_edges[ f_01_02_03 ][ 2 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 0 ] = fp_e_02_03;
      
      data_edge_firstparent_face[ e_02_23 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 1 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 1 ] = fp_e_02_23;
      
      data_edge_firstparent_face[ e_03_23 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 1 ] = fp_e_03_23;
      
      
      data_edge_firstparent_face[ e_12_13 ] = f_01_12_13;
      data_face_nextparents_of_edges[ f_01_12_13 ][ 2 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 2 ] = fp_e_12_13;
      
      data_edge_firstparent_face[ e_12_23 ] = f_02_12_23;
      data_face_nextparents_of_edges[ f_02_12_23 ][ 2 ] = fp_e_12_23;
      
      data_edge_firstparent_face[ e_13_23 ] = f_03_13_23;
      data_face_nextparents_of_edges[ f_03_13_23 ][ 2 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 2 ] = fp_e_13_23;
      
      
      data_edge_firstparent_face[ e_02_13 ] = f_01_02_13;
      data_face_nextparents_of_edges[ f_01_02_13 ][ 2 ] = f_02_03_13;
      data_face_nextparents_of_edges[ f_02_03_13 ][ 1 ] = f_02_12_13;
      data_face_nextparents_of_edges[ f_02_12_13 ][ 1 ] = f_02_13_23;
      data_face_nextparents_of_edges[ f_02_13_23 ][ 0 ] = fp_e_02_13;
      
    }
    
    
    
    
    
    
    
    /*** SET THE VERTICES AND EDGES OF ALL FACES ****/
    
    /* for each new outer face, set the new vertices */
    
    for( int f = 0; f < counter_faces; f++ )
    {
      int v00 = data_face_vertices[f][0];
      int v11 = data_face_vertices[f][1];
      int v22 = data_face_vertices[f][2];
      
      int v01 = counter_vertices + data_face_edges[f][0];
      int v02 = counter_vertices + data_face_edges[f][1];
      int v12 = counter_vertices + data_face_edges[f][2];
      
      int f_00_01_12 = 0 * counter_faces + f;
      int f_01_11_12 = 1 * counter_faces + f;
      int f_02_12_22 = 2 * counter_faces + f;
      int f_01_02_12 = 3 * counter_faces + f;
      
        
      
      
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_vertices[ f_00_01_12 ][0] = v00;
      data_face_vertices[ f_00_01_12 ][1] = v01;
      data_face_vertices[ f_00_01_12 ][2] = v02;
      
      data_face_vertices[ f_01_11_12 ][0] = v01;
      data_face_vertices[ f_01_11_12 ][1] = v11;
      data_face_vertices[ f_01_11_12 ][2] = v12;
      
      data_face_vertices[ f_02_12_22 ][0] = v02;
      data_face_vertices[ f_02_12_22 ][1] = v12;
      data_face_vertices[ f_02_12_22 ][2] = v22;
      
      data_face_vertices[ f_01_02_12 ][0] = v01;
      data_face_vertices[ f_01_02_12 ][1] = v02;
      data_face_vertices[ f_01_02_12 ][2] = v12;
      
    }
    
    /* for each new interior face, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      assert( counter_vertices <= v01 && v01 < counter_vertices + counter_edges );
      assert( counter_vertices <= v02 && v02 < counter_vertices + counter_edges );
      assert( counter_vertices <= v03 && v03 < counter_vertices + counter_edges );
      assert( counter_vertices <= v12 && v12 < counter_vertices + counter_edges );
      assert( counter_vertices <= v13 && v13 < counter_vertices + counter_edges );
      assert( counter_vertices <= v23 && v23 < counter_vertices + counter_edges );
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      data_face_vertices[ f_01_02_03 ][0] = v01;
      data_face_vertices[ f_01_02_03 ][1] = v02;
      data_face_vertices[ f_01_02_03 ][2] = v03;
      
      data_face_vertices[ f_01_12_13 ][0] = v01;
      data_face_vertices[ f_01_12_13 ][1] = v12;
      data_face_vertices[ f_01_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_23 ][0] = v02;
      data_face_vertices[ f_02_12_23 ][1] = v12;
      data_face_vertices[ f_02_12_23 ][2] = v23;
      
      data_face_vertices[ f_03_13_23 ][0] = v03;
      data_face_vertices[ f_03_13_23 ][1] = v13;
      data_face_vertices[ f_03_13_23 ][2] = v23;
      
      data_face_vertices[ f_01_02_13 ][0] = v01;
      data_face_vertices[ f_01_02_13 ][1] = v02;
      data_face_vertices[ f_01_02_13 ][2] = v13;
      
      data_face_vertices[ f_02_03_13 ][0] = v02;
      data_face_vertices[ f_02_03_13 ][1] = v03;
      data_face_vertices[ f_02_03_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_13 ][0] = v02;
      data_face_vertices[ f_02_12_13 ][1] = v12;
      data_face_vertices[ f_02_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_13_23 ][0] = v02;
      data_face_vertices[ f_02_13_23 ][1] = v13;
      data_face_vertices[ f_02_13_23 ][2] = v23;
      
      
      
    }
    
    /* for each new outer face, set the new edges */
    
    for( int f = 0; f < counter_faces; f++ )
    {
    
      int e_00_01 = 0 * counter_edges + data_face_edges[f][0];
      int e_00_02 = 0 * counter_edges + data_face_edges[f][1];
      int e_01_11 = 1 * counter_edges + data_face_edges[f][0];
      int e_02_22 = 1 * counter_edges + data_face_edges[f][1];
      int e_11_12 = 0 * counter_edges + data_face_edges[f][2];
      int e_12_22 = 1 * counter_edges + data_face_edges[f][2];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_edges[ 0 * counter_faces + f ][0] = e_00_01;
      data_face_edges[ 0 * counter_faces + f ][1] = e_00_02;
      data_face_edges[ 0 * counter_faces + f ][2] = e_01_02;
      
      data_face_edges[ 1 * counter_faces + f ][0] = e_01_11;
      data_face_edges[ 1 * counter_faces + f ][1] = e_01_12;
      data_face_edges[ 1 * counter_faces + f ][2] = e_11_12;
      
      data_face_edges[ 2 * counter_faces + f ][0] = e_02_12;
      data_face_edges[ 2 * counter_faces + f ][1] = e_02_22;
      data_face_edges[ 2 * counter_faces + f ][2] = e_12_22;
      
      data_face_edges[ 3 * counter_faces + f ][0] = e_01_02;
      data_face_edges[ 3 * counter_faces + f ][1] = e_01_12;
      data_face_edges[ 3 * counter_faces + f ][2] = e_02_12;
      
    }
    
    
    /* for each new interior face, set the new edges */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      int f_012 = data_tetrahedron_faces[t][0];
      int f_013 = data_tetrahedron_faces[t][1];
      int f_023 = data_tetrahedron_faces[t][2];
      int f_123 = data_tetrahedron_faces[t][3];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
      
      int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
      int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
      int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
      
      int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
      int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
      int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
      
      int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
      int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
      int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
      
      int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
      
      data_face_edges[ f_01_02_03 ][0] = e_01_02;
      data_face_edges[ f_01_02_03 ][1] = e_01_03;
      data_face_edges[ f_01_02_03 ][2] = e_02_03;
      
      data_face_edges[ f_01_12_13 ][0] = e_01_12;
      data_face_edges[ f_01_12_13 ][1] = e_01_13;
      data_face_edges[ f_01_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_12_23 ][0] = e_02_12;
      data_face_edges[ f_02_12_23 ][1] = e_02_23;
      data_face_edges[ f_02_12_23 ][2] = e_12_23;
      
      data_face_edges[ f_03_13_23 ][0] = e_03_13;
      data_face_edges[ f_03_13_23 ][1] = e_03_23;
      data_face_edges[ f_03_13_23 ][2] = e_13_23;
      
      data_face_edges[ f_01_02_13 ][0] = e_01_02;
      data_face_edges[ f_01_02_13 ][1] = e_01_13;
      data_face_edges[ f_01_02_13 ][2] = e_02_13;
      
      data_face_edges[ f_02_03_13 ][0] = e_02_03;
      data_face_edges[ f_02_03_13 ][1] = e_02_13;
      data_face_edges[ f_02_03_13 ][2] = e_03_13;
      
      data_face_edges[ f_02_12_13 ][0] = e_02_12;
      data_face_edges[ f_02_12_13 ][1] = e_02_13;
      data_face_edges[ f_02_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_13_23 ][0] = e_02_13;
      data_face_edges[ f_02_13_23 ][1] = e_02_23;
      data_face_edges[ f_02_13_23 ][2] = e_13_23;
      
    }
    
    
    
    /**************************/
    /***   ALL DIMENSIONS  ****/
    /**************************/
    
    
    /*** TREAT THE OLD VERTICES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    /* 1. for each old vertex, set the new parent tetrahedron */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      
      int p = data_vertex_firstparent_tetrahedron[v];
      
      assert( p != nullindex && 0 <= p && p < counter_tetrahedra );
      
      int vi = nullindex;
      if( data_tetrahedron_vertices[p][0] == v ) vi = 0;
      if( data_tetrahedron_vertices[p][1] == v ) vi = 1;
      if( data_tetrahedron_vertices[p][2] == v ) vi = 2;
      if( data_tetrahedron_vertices[p][3] == v ) vi = 3;
      assert( vi != nullindex );
      
      assert( data_tetrahedron_vertices[p][vi] == v );
      
      data_vertex_firstparent_tetrahedron[v] = vi * counter_tetrahedra + p;
      
    }
    
    
    /* 2. for each old tetrahedron, relocate the data of the old vertices' parent tetrahedron */
    
    for( int t  = 0; t  < counter_tetrahedra;  t++ )
    for( int vi = 0; vi <                  4; vi++ )
    {
      
      int q = data_tetrahedron_nextparents_of_vertices[t][vi];
      
      int v = data_tetrahedron_vertices[t][vi];
      
      if( q == nullindex ) {
        
        data_tetrahedron_nextparents_of_vertices[ vi * counter_tetrahedra + t ][vi] = nullindex;
        
      } else {
        
        int vinp = nullindex;
        if( data_tetrahedron_vertices[q][0] == v ) vinp = 0;
        if( data_tetrahedron_vertices[q][1] == v ) vinp = 1;
        if( data_tetrahedron_vertices[q][2] == v ) vinp = 2;
        if( data_tetrahedron_vertices[q][3] == v ) vinp = 3;
        assert( vinp != nullindex );
        
        assert( data_tetrahedron_vertices[q][vinp] == v );
        
        data_tetrahedron_nextparents_of_vertices[ vi * counter_tetrahedra + t ][vi] = vinp * counter_tetrahedra + q;
      
      } 
      
    }
    
    /*
    *       00 01 02 03
    *       01 11 12 13
    *       02 12 22 23
    *       03 13 23 33
    *       
    *       01 02 03 13
    *       01 02 12 13 
    *       02 03 13 23
    *       02 12 13 23
    */
       
        
    /*** TREAT THE NEW VERTICES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int v01 = counter_vertices + data_tetrahedron_edges[ t ][ 0 ];
        int v02 = counter_vertices + data_tetrahedron_edges[ t ][ 1 ];
        int v03 = counter_vertices + data_tetrahedron_edges[ t ][ 2 ];
        int v12 = counter_vertices + data_tetrahedron_edges[ t ][ 3 ];
        int v13 = counter_vertices + data_tetrahedron_edges[ t ][ 4 ];
        int v23 = counter_vertices + data_tetrahedron_edges[ t ][ 5 ];
        
        int fp_v01 = data_vertex_firstparent_tetrahedron[ v01 ];
        int fp_v02 = data_vertex_firstparent_tetrahedron[ v02 ];
        int fp_v03 = data_vertex_firstparent_tetrahedron[ v03 ];
        int fp_v12 = data_vertex_firstparent_tetrahedron[ v12 ];
        int fp_v13 = data_vertex_firstparent_tetrahedron[ v13 ];
        int fp_v23 = data_vertex_firstparent_tetrahedron[ v23 ];
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        data_vertex_firstparent_tetrahedron[ v01 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 1 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 0 ] = fp_v01;
        
        data_vertex_firstparent_tetrahedron[ v02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 2 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 1 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 1 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 0 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 0 ] = fp_v02;
        
        data_vertex_firstparent_tetrahedron[ v03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_vertices[ t_00_01_02_03 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 0 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 1 ] = fp_v03;
        
        data_vertex_firstparent_tetrahedron[ v12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 2 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 1 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 1 ] = fp_v12;
        
        data_vertex_firstparent_tetrahedron[ v13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_11_12_13 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 1 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_03_13 ][ 3 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_vertices[ t_01_02_12_13 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 2 ] = fp_v13;
        
        data_vertex_firstparent_tetrahedron[ v23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_22_23 ][ 3 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_vertices[ t_03_13_23_33 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_03_13_23 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_vertices[ t_02_12_13_23 ][ 3 ] = fp_v23;
        
    }
    
    /*** TREAT THE BISECTED EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int e_00_01 = data_tetrahedron_edges[ t ][ 0 ];
        int e_00_02 = data_tetrahedron_edges[ t ][ 1 ];
        int e_00_03 = data_tetrahedron_edges[ t ][ 2 ];
        int e_11_12 = data_tetrahedron_edges[ t ][ 3 ];
        int e_11_13 = data_tetrahedron_edges[ t ][ 4 ];
        int e_22_23 = data_tetrahedron_edges[ t ][ 5 ];
        
        int e_01_11 = counter_edges + data_tetrahedron_edges[ t ][ 0 ];
        int e_02_22 = counter_edges + data_tetrahedron_edges[ t ][ 1 ];
        int e_03_33 = counter_edges + data_tetrahedron_edges[ t ][ 2 ];
        int e_12_22 = counter_edges + data_tetrahedron_edges[ t ][ 3 ];
        int e_13_33 = counter_edges + data_tetrahedron_edges[ t ][ 4 ];
        int e_23_33 = counter_edges + data_tetrahedron_edges[ t ][ 5 ];
        
        
        int fp_e_00_01 = data_edge_firstparent_tetrahedron[ e_00_01 ];
        int fp_e_00_02 = data_edge_firstparent_tetrahedron[ e_00_02 ];
        int fp_e_00_03 = data_edge_firstparent_tetrahedron[ e_00_03 ];
        int fp_e_11_12 = data_edge_firstparent_tetrahedron[ e_11_12 ];
        int fp_e_11_13 = data_edge_firstparent_tetrahedron[ e_11_13 ];
        int fp_e_22_23 = data_edge_firstparent_tetrahedron[ e_22_23 ];
        
        int fp_e_01_11 = data_edge_firstparent_tetrahedron[ e_01_11 ];
        int fp_e_02_22 = data_edge_firstparent_tetrahedron[ e_02_22 ];
        int fp_e_03_33 = data_edge_firstparent_tetrahedron[ e_03_33 ];
        int fp_e_12_22 = data_edge_firstparent_tetrahedron[ e_12_22 ];
        int fp_e_13_33 = data_edge_firstparent_tetrahedron[ e_13_33 ];
        int fp_e_23_33 = data_edge_firstparent_tetrahedron[ e_23_33 ];
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;

        data_edge_firstparent_tetrahedron[ e_00_01 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 0 ] = fp_e_00_01;
        
        data_edge_firstparent_tetrahedron[ e_00_02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 1 ] = fp_e_00_02;
        
        data_edge_firstparent_tetrahedron[ e_00_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 2 ] = fp_e_00_03;
        
        data_edge_firstparent_tetrahedron[ e_11_12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 3 ] = fp_e_11_12;
        
        data_edge_firstparent_tetrahedron[ e_11_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 4 ] = fp_e_11_13;
        
        data_edge_firstparent_tetrahedron[ e_22_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 5 ] = fp_e_22_23;
        
        
        data_edge_firstparent_tetrahedron[ e_01_11 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 0 ] = fp_e_01_11;
        
        data_edge_firstparent_tetrahedron[ e_02_22 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 1 ] = fp_e_02_22;
        
        data_edge_firstparent_tetrahedron[ e_03_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 2 ] = fp_e_03_33;
        
        data_edge_firstparent_tetrahedron[ e_12_22 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 3 ] = fp_e_12_22;
        
        data_edge_firstparent_tetrahedron[ e_13_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 4 ] = fp_e_13_33;
        
        data_edge_firstparent_tetrahedron[ e_23_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 5 ] = fp_e_23_33;
                
        
    }
    
    
    
    
    /*** TREAT THE FACE-BASED EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    /*** TREAT THE SINGLE INTERIOR EDGES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int f_012 = data_tetrahedron_faces[ t ][ 0 ];
        int f_013 = data_tetrahedron_faces[ t ][ 1 ];
        int f_023 = data_tetrahedron_faces[ t ][ 2 ];
        int f_123 = data_tetrahedron_faces[ t ][ 3 ];
        
        /* 01 02 12 */
        /* 01 03 13 */
        /* 02 03 23 */
        /* 12 13 23 */
        
        int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
        int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
        int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
        
        int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
        int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
        int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
        
        int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
        int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
        int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
        
        int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
        int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
        int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
        
        int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
    
        int fp_e_01_02 = data_edge_firstparent_tetrahedron[ e_01_02 ];
        int fp_e_01_12 = data_edge_firstparent_tetrahedron[ e_01_12 ];
        int fp_e_02_12 = data_edge_firstparent_tetrahedron[ e_02_12 ];
        
        int fp_e_01_03 = data_edge_firstparent_tetrahedron[ e_01_03 ];
        int fp_e_01_13 = data_edge_firstparent_tetrahedron[ e_01_13 ];
        int fp_e_03_13 = data_edge_firstparent_tetrahedron[ e_03_13 ];
        
        int fp_e_02_03 = data_edge_firstparent_tetrahedron[ e_02_03 ];
        int fp_e_02_23 = data_edge_firstparent_tetrahedron[ e_02_23 ];
        int fp_e_03_23 = data_edge_firstparent_tetrahedron[ e_03_23 ];
        
        int fp_e_12_13 = data_edge_firstparent_tetrahedron[ e_12_13 ];
        int fp_e_12_23 = data_edge_firstparent_tetrahedron[ e_12_23 ];
        int fp_e_13_23 = data_edge_firstparent_tetrahedron[ e_13_23 ];
        
        int fp_e_02_13 = data_edge_firstparent_tetrahedron[ e_02_13 ];
        assert( fp_e_02_13 == nullindex );
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        
        data_edge_firstparent_tetrahedron[ e_01_02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 3 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 0 ] = fp_e_01_02;
        
        data_edge_firstparent_tetrahedron[ e_01_12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 2 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 2 ] = fp_e_01_12;
        
        data_edge_firstparent_tetrahedron[ e_02_12 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 0 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 0 ] = fp_e_02_12;
        
        
        data_edge_firstparent_tetrahedron[ e_01_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 4 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 1 ] = fp_e_01_03;
        
        data_edge_firstparent_tetrahedron[ e_01_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 2 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 2 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 2 ] = fp_e_01_13;
        
        data_edge_firstparent_tetrahedron[ e_03_13 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 0 ] = fp_e_03_13;
        
        
        data_edge_firstparent_tetrahedron[ e_02_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_edges[ t_00_01_02_03 ][ 5 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 0 ] = fp_e_02_03;
        
        data_edge_firstparent_tetrahedron[ e_02_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 2 ] = fp_e_02_23;
        
        data_edge_firstparent_tetrahedron[ e_03_23 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 2 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 3 ] = fp_e_03_23;
        
        
        data_edge_firstparent_tetrahedron[ e_12_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_11_12_13 ][ 5 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 5 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 3 ] = fp_e_12_13;
        
        data_edge_firstparent_tetrahedron[ e_12_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_22_23 ][ 4 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 4 ] = fp_e_12_23;
        
        data_edge_firstparent_tetrahedron[ e_13_23 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_edges[ t_03_13_23_33 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 5 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 5 ] = fp_e_13_23;
        
        
        
        data_edge_firstparent_tetrahedron[ e_02_13 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_03_13 ][ 4 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_edges[ t_01_02_12_13 ][ 4 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_03_13_23 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_edges[ t_02_12_13_23 ][ 1 ] = fp_e_02_13;
        
        
        
        
    }
    
    
    
    
    
    /*** TREAT THE OUTER FACES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    /* TODO */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int f_012 = data_tetrahedron_faces[ t ][ 0 ];
        int f_013 = data_tetrahedron_faces[ t ][ 1 ];
        int f_023 = data_tetrahedron_faces[ t ][ 2 ];
        int f_123 = data_tetrahedron_faces[ t ][ 3 ];
        
        
        int f_00_01_02 = 0 * counter_faces + f_012;
        int f_01_11_12 = 1 * counter_faces + f_012;
        int f_02_12_22 = 2 * counter_faces + f_012;
        int f_01_02_12 = 3 * counter_faces + f_012;
        
        int f_00_01_03 = 0 * counter_faces + f_013;
        int f_01_11_13 = 1 * counter_faces + f_013;
        int f_03_13_33 = 2 * counter_faces + f_013;
        int f_01_03_13 = 3 * counter_faces + f_013;
        
        int f_00_02_03 = 0 * counter_faces + f_023;
        int f_02_22_23 = 1 * counter_faces + f_023;
        int f_03_23_33 = 2 * counter_faces + f_023;
        int f_02_03_23 = 3 * counter_faces + f_023;
        
        int f_11_12_13 = 0 * counter_faces + f_123;
        int f_12_22_23 = 1 * counter_faces + f_123;
        int f_13_23_33 = 2 * counter_faces + f_123;
        int f_12_13_23 = 3 * counter_faces + f_123;
        
        
        int fp_f_00_01_02 = data_face_firstparent_tetrahedron[ f_00_01_02 ];
        int fp_f_01_11_12 = data_face_firstparent_tetrahedron[ f_01_11_12 ];
        int fp_f_02_12_22 = data_face_firstparent_tetrahedron[ f_02_12_22 ];
        int fp_f_01_02_12 = data_face_firstparent_tetrahedron[ f_01_02_12 ];
        
        int fp_f_00_01_03 = data_face_firstparent_tetrahedron[ f_00_01_03 ];
        int fp_f_01_11_13 = data_face_firstparent_tetrahedron[ f_01_11_13 ];
        int fp_f_03_13_33 = data_face_firstparent_tetrahedron[ f_03_13_33 ];
        int fp_f_01_03_13 = data_face_firstparent_tetrahedron[ f_01_03_13 ];
        
        int fp_f_00_02_03 = data_face_firstparent_tetrahedron[ f_00_02_03 ];
        int fp_f_02_22_23 = data_face_firstparent_tetrahedron[ f_02_22_23 ];
        int fp_f_03_23_33 = data_face_firstparent_tetrahedron[ f_03_23_33 ];
        int fp_f_02_03_23 = data_face_firstparent_tetrahedron[ f_02_03_23 ];
        
        int fp_f_11_12_13 = data_face_firstparent_tetrahedron[ f_11_12_13 ];
        int fp_f_12_22_23 = data_face_firstparent_tetrahedron[ f_12_22_23 ];
        int fp_f_13_23_33 = data_face_firstparent_tetrahedron[ f_13_23_33 ];
        int fp_f_12_13_23 = data_face_firstparent_tetrahedron[ f_12_13_23 ];
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        data_face_firstparent_tetrahedron[ f_00_01_02 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_faces[ t_00_01_02_03 ][ 0 ] = fp_f_00_01_02;
        
        data_face_firstparent_tetrahedron[ f_01_11_12 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_11_12_13 ][ 0 ] = fp_f_01_11_12;
        
        data_face_firstparent_tetrahedron[ f_02_12_22 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_22_23 ][ 0 ] = fp_f_02_12_22;
        
        data_face_firstparent_tetrahedron[ f_01_02_12 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_12_13 ][ 0 ] = fp_f_01_02_12;
        
        
        data_face_firstparent_tetrahedron[ f_00_01_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_faces[ t_00_01_02_03 ][ 1 ] = fp_f_00_01_03;
        
        data_face_firstparent_tetrahedron[ f_01_11_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_11_12_13 ][ 1 ] = fp_f_01_11_13;
        
        data_face_firstparent_tetrahedron[ f_03_13_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_faces[ t_03_13_23_33 ][ 1 ] = fp_f_03_13_33;
        
        data_face_firstparent_tetrahedron[ f_01_03_13 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_03_13 ][ 2 ] = fp_f_01_03_13;
        
        
        data_face_firstparent_tetrahedron[ f_00_02_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_faces[ t_00_01_02_03 ][ 2 ] = fp_f_00_02_03;
        
        data_face_firstparent_tetrahedron[ f_02_22_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_22_23 ][ 2 ] = fp_f_02_22_23;
        
        data_face_firstparent_tetrahedron[ f_03_23_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_faces[ t_03_13_23_33 ][ 2 ] = fp_f_03_23_33;
        
        data_face_firstparent_tetrahedron[ f_02_03_23 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_03_13_23 ][ 1 ] = fp_f_02_03_23;
        
        
        data_face_firstparent_tetrahedron[ f_11_12_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_11_12_13 ][ 3 ] = fp_f_11_12_13;
        
        data_face_firstparent_tetrahedron[ f_12_22_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_22_23 ][ 3 ] = fp_f_12_22_23;
        
        data_face_firstparent_tetrahedron[ f_13_23_33 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_faces[ t_03_13_23_33 ][ 3 ] = fp_f_13_23_33;
        
        data_face_firstparent_tetrahedron[ f_12_13_23 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_13_23 ][ 3 ] = fp_f_12_13_23;
        
    }
    
    
    /*** TREAT THE INNER FACES AND THEIR CONNECTION TO TETRAHEDRA ****/
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        
        int f_012 = data_tetrahedron_faces[ t ][ 0 ];
        int f_013 = data_tetrahedron_faces[ t ][ 1 ];
        int f_023 = data_tetrahedron_faces[ t ][ 2 ];
        int f_123 = data_tetrahedron_faces[ t ][ 3 ];
        
        
        int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
        int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
        int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
        int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
        int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
        int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
        int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
        int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
        
        int fp_f_01_02_03 = data_face_firstparent_tetrahedron[ f_01_02_03 ];
        int fp_f_01_12_13 = data_face_firstparent_tetrahedron[ f_01_12_13 ];
        int fp_f_02_12_23 = data_face_firstparent_tetrahedron[ f_02_12_23 ];
        int fp_f_03_13_23 = data_face_firstparent_tetrahedron[ f_03_13_23 ];
        int fp_f_01_02_13 = data_face_firstparent_tetrahedron[ f_01_02_13 ];
        int fp_f_02_03_13 = data_face_firstparent_tetrahedron[ f_02_03_13 ];
        int fp_f_02_12_13 = data_face_firstparent_tetrahedron[ f_02_12_13 ];
        int fp_f_02_13_23 = data_face_firstparent_tetrahedron[ f_02_13_23 ];
        
        /*
        * 01 02 03
        * 01 12 13
        * 02 12 23
        * 03 13 23
        *       
        * 01 02 13
        * 02 03 13
        * 02 12 13
        * 02 13 23
        */
        
        
        int t_00_01_02_03 = 0 * counter_tetrahedra + t;
        int t_01_11_12_13 = 1 * counter_tetrahedra + t;
        int t_02_12_22_23 = 2 * counter_tetrahedra + t;
        int t_03_13_23_33 = 3 * counter_tetrahedra + t;
        int t_01_02_03_13 = 4 * counter_tetrahedra + t;
        int t_01_02_12_13 = 5 * counter_tetrahedra + t;
        int t_02_03_13_23 = 6 * counter_tetrahedra + t;
        int t_02_12_13_23 = 7 * counter_tetrahedra + t;
        
        data_face_firstparent_tetrahedron[ f_01_02_03 ] = t_00_01_02_03;
        data_tetrahedron_nextparents_of_faces[ t_00_01_02_03 ][ 3 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_03_13 ][ 0 ] = fp_f_01_02_03;
        
        data_face_firstparent_tetrahedron[ f_01_12_13 ] = t_01_11_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_11_12_13 ][ 2 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_12_13 ][ 2 ] = fp_f_01_12_13;
        
        data_face_firstparent_tetrahedron[ f_02_12_23 ] = t_02_12_22_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_22_23 ][ 1 ] = fp_f_02_12_23;
        
        data_face_firstparent_tetrahedron[ f_03_13_23 ] = t_03_13_23_33;
        data_tetrahedron_nextparents_of_faces[ t_03_13_23_33 ][ 0 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_03_13_23 ][ 3 ] = fp_f_03_13_23;
        
        data_face_firstparent_tetrahedron[ f_01_02_13 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_03_13 ][ 1 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_12_13 ][ 1 ] = fp_f_01_02_13;
        
        data_face_firstparent_tetrahedron[ f_02_03_13 ] = t_01_02_03_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_03_13 ][ 3 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_03_13_23 ][ 2 ] = fp_f_02_03_13;
        
        data_face_firstparent_tetrahedron[ f_02_12_13 ] = t_01_02_12_13;
        data_tetrahedron_nextparents_of_faces[ t_01_02_12_13 ][ 3 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_13_23 ][ 0 ] = fp_f_02_12_13;
        
        data_face_firstparent_tetrahedron[ f_02_13_23 ] = t_02_03_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_03_13_23 ][ 2 ] = t_02_12_13_23;
        data_tetrahedron_nextparents_of_faces[ t_02_12_13_23 ][ 2 ] = fp_f_02_13_23;
        
        
        
    }
    
    
    
    
    
    
    
    
    /* for each new tetrahedron, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v00 = data_tetrahedron_vertices[t][0];
      int v11 = data_tetrahedron_vertices[t][1];
      int v22 = data_tetrahedron_vertices[t][2];
      int v33 = data_tetrahedron_vertices[t][3];
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      //       00 01 02 03
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][0] = v00;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][1] = v01;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][2] = v02;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][3] = v03;
      
      //       01 11 12 13
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][1] = v11;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][3] = v13;
      
      //       02 12 22 23
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][2] = v22;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][3] = v23;
      
      //       03 13 23 33
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][0] = v03;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][1] = v13;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][2] = v23;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][3] = v33;
      
      
      //       01 02 03 13
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][2] = v03;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][3] = v13;
      
      //       01 02 12 13 
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][3] = v13;
      
      //       02 03 13 23
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][1] = v03;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][3] = v23;
      
      //       02 12 13 23
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][3] = v23;
      
    }
    
    /* for each new tetrahedron, set the new edges */
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        // 00 11 22 33 
        
        // From old edges ...
        int e_00_01 = data_tetrahedron_edges[t][0] + 0 * counter_edges;
        int e_01_11 = data_tetrahedron_edges[t][0] + 1 * counter_edges;
        
        int e_00_02 = data_tetrahedron_edges[t][1] + 0 * counter_edges;
        int e_02_22 = data_tetrahedron_edges[t][1] + 1 * counter_edges;
        
        int e_00_03 = data_tetrahedron_edges[t][2] + 0 * counter_edges;
        int e_03_33 = data_tetrahedron_edges[t][2] + 1 * counter_edges;
        
        int e_11_12 = data_tetrahedron_edges[t][3] + 0 * counter_edges;
        int e_12_22 = data_tetrahedron_edges[t][3] + 1 * counter_edges;
        
        int e_11_13 = data_tetrahedron_edges[t][4] + 0 * counter_edges;
        int e_13_33 = data_tetrahedron_edges[t][4] + 1 * counter_edges;
        
        int e_22_23 = data_tetrahedron_edges[t][5] + 0 * counter_edges;
        int e_23_33 = data_tetrahedron_edges[t][5] + 1 * counter_edges;
        
        
        // ... from old faces ... 
        
        // 00 11 22  
        int e_01_02 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][0];
        int e_01_12 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][0];
        int e_02_12 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][0];
        
        // 00 11 33 
        int e_01_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][1];
        int e_01_13 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][1];
        int e_03_13 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][1];
        
        // 00 22 33 
        int e_02_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][2];
        int e_02_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][2];
        int e_03_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][2];
        
        // 11 22 33 
        int e_12_13 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][3];
        int e_12_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][3];
        int e_13_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][3];
        
        // ... and the single new internal one. 
        
        int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
        
        
        // fill in to the tetrahedra
        
        //       00 01 02 03
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][0] = e_00_01;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][1] = e_00_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][2] = e_00_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][3] = e_01_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][4] = e_01_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][5] = e_02_03;
        
        //       01 11 12 13
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][0] = e_01_11;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][1] = e_01_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][3] = e_11_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][4] = e_11_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][5] = e_12_13;
        
        //       02 12 22 23
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][0] = e_02_12;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][1] = e_02_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][3] = e_12_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][4] = e_12_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][5] = e_22_23;
        
        //       03 13 23 33
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][0] = e_03_13;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][1] = e_03_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][2] = e_03_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][3] = e_13_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][4] = e_13_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][5] = e_23_33;
        
        
        //       01 02 03 13
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][0] = e_01_02;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][1] = e_01_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][3] = e_02_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][4] = e_02_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][5] = e_03_13;
        
        //       01 02 12 13 
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][0] = e_01_02;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][1] = e_01_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][3] = e_02_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][4] = e_02_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][5] = e_12_13;
        
        //       02 03 13 23
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][0] = e_02_03;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][1] = e_02_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][3] = e_03_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][4] = e_03_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][5] = e_13_23;
        
        //       02 12 13 23
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][0] = e_02_12;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][1] = e_02_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][3] = e_12_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][4] = e_12_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][5] = e_13_23;
        
    }
    
    /* for each new tetrahedron, set the new faces */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      // 00 11 22 33 
      
      // From old faces ... 
      
      // 00 11 22 
      int f_00_01_02 = data_tetrahedron_faces[t][0] + 0 * counter_faces;
      int f_01_11_12 = data_tetrahedron_faces[t][0] + 1 * counter_faces;
      int f_02_12_22 = data_tetrahedron_faces[t][0] + 2 * counter_faces;
      int f_01_02_12 = data_tetrahedron_faces[t][0] + 3 * counter_faces;
      
      // 00 11 33 
      int f_00_01_03 = data_tetrahedron_faces[t][1] + 0 * counter_faces;
      int f_01_11_13 = data_tetrahedron_faces[t][1] + 1 * counter_faces;
      int f_03_13_33 = data_tetrahedron_faces[t][1] + 2 * counter_faces;
      int f_01_03_13 = data_tetrahedron_faces[t][1] + 3 * counter_faces;
      
      // 00 22 33 
      int f_00_02_03 = data_tetrahedron_faces[t][2] + 0 * counter_faces;
      int f_02_22_23 = data_tetrahedron_faces[t][2] + 1 * counter_faces;
      int f_03_23_33 = data_tetrahedron_faces[t][2] + 2 * counter_faces;
      int f_02_03_23 = data_tetrahedron_faces[t][2] + 3 * counter_faces;
      
      // 11 22 33 
      int f_11_12_13 = data_tetrahedron_faces[t][3] + 0 * counter_faces;
      int f_12_22_23 = data_tetrahedron_faces[t][3] + 1 * counter_faces;
      int f_13_23_33 = data_tetrahedron_faces[t][3] + 2 * counter_faces;
      int f_12_13_23 = data_tetrahedron_faces[t][3] + 3 * counter_faces;
      
      // ... new internal faces. 
      
      // 01 02 03
      // 01 12 13
      // 02 12 23
      // 03 13 23
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      
      // 01 02 13
      // 02 03 13
      // 02 12 13
      // 02 13 23
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      
      
      //       00 01 02 03
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][0] = f_00_01_02;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][1] = f_00_01_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][2] = f_00_02_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][3] = f_01_02_03;
      
      //       01 11 12 13
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][0] = f_01_11_12;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][1] = f_01_11_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][2] = f_01_12_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][3] = f_11_12_13;
      
      //       02 12 22 23
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][0] = f_02_12_22;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][1] = f_02_12_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][2] = f_02_22_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][3] = f_12_22_23;
      
      //       03 13 23 33
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][0] = f_03_13_23;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][1] = f_03_13_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][2] = f_03_23_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][3] = f_13_23_33;
      
      
      //       01 02 03 13
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][0] = f_01_02_03;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][1] = f_01_02_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][2] = f_01_03_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][3] = f_02_03_13;
      
      //       01 02 12 13 
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][0] = f_01_02_12;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][1] = f_01_02_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][2] = f_01_12_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][3] = f_02_12_13;
      
      //       02 03 13 23
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][0] = f_02_03_13;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][1] = f_02_03_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][2] = f_02_13_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][3] = f_03_13_23;
      
      //       02 12 13 23
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][0] = f_02_12_13;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][1] = f_02_12_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][2] = f_02_13_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][3] = f_12_13_23;
      
    }
    
    
    
    
    /* update the counters */
    
    counter_vertices   = new_counter_vertices;
    counter_edges      = new_counter_edges;
    counter_faces      = new_counter_faces;
    counter_tetrahedra = new_counter_tetrahedra;
    
    /* DONE */
    
//     check();
}








