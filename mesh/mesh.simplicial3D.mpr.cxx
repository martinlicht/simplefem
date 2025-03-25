

#include "../base/include.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"

#include "mesh.simplicial3D.hpp"










void MeshSimplicial3D::midpoint_refinement( const int t )
{
    // check();
    
    assert( 0 <= t && t < counter_tetrahedra );
    
    /* Allocate memory */
    
    data_tetrahedron_faces.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    data_face_firstparent_tetrahedron.resize    ( counter_faces + 6,                                         nullindex   );
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    
    data_tetrahedron_edges.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_edge_firstparent_tetrahedron.resize    ( counter_edges + 4,                                                               nullindex   );
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    data_tetrahedron_vertices.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                      nullindex   );
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_edges.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
    data_edge_firstparent_face.resize    ( counter_edges + 4,                         nullindex   );
    data_face_nextparents_of_edges.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
    
    data_face_vertices.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_face.resize     ( counter_vertices + 1,                      nullindex   );
    data_face_nextparents_of_vertices.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
    
    data_edge_vertices.resize               ( counter_edges + 4,   { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,             nullindex   );
    data_edge_nextparents_of_vertices.resize( counter_edges + 4,   { nullindex, nullindex } );

    auto t_flag = flags_tetrahedra[t];
    flags_tetrahedra.resize( counter_tetrahedra + 3, t_flag );
    flags_faces.resize     ( counter_faces + 6,      t_flag );
    flags_edges.resize     ( counter_edges + 4,      t_flag );
    flags_vertices.resize  ( counter_vertices + 1,   t_flag );

    
    getCoordinates().addcoordinates( 1 );
    
    
    // 0. Check the array sizes 
    // 1. check that the internal data of each simplex make sense on each dimension 
    // 2. check the uniqueness of simplices 
    // 3. check the data of each simplex accross dimensions 
    
    /* load the new coordinate */
    
    getCoordinates().setdata_by_vertex( counter_vertices, get_tetrahedron_midpoint( t ) );
    
    /* assemble the data and auxiliary variables */

    // old subsimplices 

    // 0 1 2, 0 1 3, 0 2 3, 1 2 3
    const int t_f0_012 = data_tetrahedron_faces[t][0];
    const int t_f1_013 = data_tetrahedron_faces[t][1];
    const int t_f2_023 = data_tetrahedron_faces[t][2];
    const int t_f3_123 = data_tetrahedron_faces[t][3];
    
    // 0 1, 0 2, 0 3, 1 2, 1 3, 2 3
    const int t_e0_01 = data_tetrahedron_edges[t][0];
    const int t_e1_02 = data_tetrahedron_edges[t][1];
    const int t_e2_03 = data_tetrahedron_edges[t][2];
    const int t_e3_12 = data_tetrahedron_edges[t][3];
    const int t_e4_13 = data_tetrahedron_edges[t][4];
    const int t_e5_23 = data_tetrahedron_edges[t][5];
    
    // 0, 1, 2, 3
    const int t_v0 = data_tetrahedron_vertices[t][0];
    const int t_v1 = data_tetrahedron_vertices[t][1];
    const int t_v2 = data_tetrahedron_vertices[t][2];
    const int t_v3 = data_tetrahedron_vertices[t][3];
    
    // old subsimplices next parents

    const int x_t_f0_012 = data_tetrahedron_nextparents_of_faces[t][0];
    const int x_t_f1_013 = data_tetrahedron_nextparents_of_faces[t][1];
    const int x_t_f2_023 = data_tetrahedron_nextparents_of_faces[t][2];
    const int x_t_f3_123 = data_tetrahedron_nextparents_of_faces[t][3];
    
    const int x_t_e0_01 = data_tetrahedron_nextparents_of_edges[t][0];
    const int x_t_e1_02 = data_tetrahedron_nextparents_of_edges[t][1];
    const int x_t_e2_03 = data_tetrahedron_nextparents_of_edges[t][2];
    const int x_t_e3_12 = data_tetrahedron_nextparents_of_edges[t][3];
    const int x_t_e4_13 = data_tetrahedron_nextparents_of_edges[t][4];
    const int x_t_e5_23 = data_tetrahedron_nextparents_of_edges[t][5];
    
    const int x_t_v0 = data_tetrahedron_nextparents_of_vertices[t][0];
    const int x_t_v1 = data_tetrahedron_nextparents_of_vertices[t][1];
    const int x_t_v2 = data_tetrahedron_nextparents_of_vertices[t][2];
    const int x_t_v3 = data_tetrahedron_nextparents_of_vertices[t][3];
    
    // old subsimplices next parents

    const int fp_e0_01_f = data_edge_firstparent_face[t_e0_01];
    const int fp_e1_02_f = data_edge_firstparent_face[t_e1_02];
    const int fp_e2_03_f = data_edge_firstparent_face[t_e2_03];
    const int fp_e3_12_f = data_edge_firstparent_face[t_e3_12];
    const int fp_e4_13_f = data_edge_firstparent_face[t_e4_13];
    const int fp_e5_23_f = data_edge_firstparent_face[t_e5_23];
    
    const int fp_v0_f = data_vertex_firstparent_face[t_v0];
    const int fp_v1_f = data_vertex_firstparent_face[t_v1];
    const int fp_v2_f = data_vertex_firstparent_face[t_v2];
    const int fp_v3_f = data_vertex_firstparent_face[t_v3];
    
    const int fp_v0_e = data_vertex_firstparent_edge[t_v0];
    const int fp_v1_e = data_vertex_firstparent_edge[t_v1];
    const int fp_v2_e = data_vertex_firstparent_edge[t_v2];
    const int fp_v3_e = data_vertex_firstparent_edge[t_v3];
    
    
    // indices of the new subsimplices

    const int t0_n012 = t;
    const int t1_n013 = counter_tetrahedra;
    const int t2_n023 = counter_tetrahedra + 1;
    const int t3_n123 = counter_tetrahedra + 2;
    
    const int f0_n01 = counter_faces + 0;
    const int f1_n02 = counter_faces + 1;
    const int f2_n03 = counter_faces + 2;
    const int f3_n12 = counter_faces + 3;
    const int f4_n13 = counter_faces + 4;
    const int f5_n23 = counter_faces + 5;
    
    const int e0_n0 = counter_edges + 0;
    const int e1_n1 = counter_edges + 1;
    const int e2_n2 = counter_edges + 2;
    const int e3_n3 = counter_edges + 3;
    
    const int vn = counter_vertices;
    
    
    /* Fill in the data how to go downwards */ 

    data_tetrahedron_vertices[ t0_n012 ][0] = vn;
    data_tetrahedron_vertices[ t0_n012 ][1] = t_v0;
    data_tetrahedron_vertices[ t0_n012 ][2] = t_v1;
    data_tetrahedron_vertices[ t0_n012 ][3] = t_v2;
    
    data_tetrahedron_vertices[ t1_n013 ][0] = vn;
    data_tetrahedron_vertices[ t1_n013 ][1] = t_v0;
    data_tetrahedron_vertices[ t1_n013 ][2] = t_v1;
    data_tetrahedron_vertices[ t1_n013 ][3] = t_v3;
    
    data_tetrahedron_vertices[ t2_n023 ][0] = vn;
    data_tetrahedron_vertices[ t2_n023 ][1] = t_v0;
    data_tetrahedron_vertices[ t2_n023 ][2] = t_v2;
    data_tetrahedron_vertices[ t2_n023 ][3] = t_v3;
    
    data_tetrahedron_vertices[ t3_n123 ][0] = vn;
    data_tetrahedron_vertices[ t3_n123 ][1] = t_v1;
    data_tetrahedron_vertices[ t3_n123 ][2] = t_v2;
    data_tetrahedron_vertices[ t3_n123 ][3] = t_v3;
    
    data_tetrahedron_edges[ t0_n012 ][0] = e0_n0;
    data_tetrahedron_edges[ t0_n012 ][1] = e1_n1;
    data_tetrahedron_edges[ t0_n012 ][2] = e2_n2;
    data_tetrahedron_edges[ t0_n012 ][3] = t_e0_01;
    data_tetrahedron_edges[ t0_n012 ][4] = t_e1_02;
    data_tetrahedron_edges[ t0_n012 ][5] = t_e3_12;
    
    data_tetrahedron_edges[ t1_n013 ][0] = e0_n0;
    data_tetrahedron_edges[ t1_n013 ][1] = e1_n1;
    data_tetrahedron_edges[ t1_n013 ][2] = e3_n3;
    data_tetrahedron_edges[ t1_n013 ][3] = t_e0_01;
    data_tetrahedron_edges[ t1_n013 ][4] = t_e2_03;
    data_tetrahedron_edges[ t1_n013 ][5] = t_e4_13;
    
    data_tetrahedron_edges[ t2_n023 ][0] = e0_n0;
    data_tetrahedron_edges[ t2_n023 ][1] = e2_n2;
    data_tetrahedron_edges[ t2_n023 ][2] = e3_n3;
    data_tetrahedron_edges[ t2_n023 ][3] = t_e1_02;
    data_tetrahedron_edges[ t2_n023 ][4] = t_e2_03;
    data_tetrahedron_edges[ t2_n023 ][5] = t_e5_23;
    
    data_tetrahedron_edges[ t3_n123 ][0] = e1_n1;
    data_tetrahedron_edges[ t3_n123 ][1] = e2_n2;
    data_tetrahedron_edges[ t3_n123 ][2] = e3_n3;
    data_tetrahedron_edges[ t3_n123 ][3] = t_e3_12;
    data_tetrahedron_edges[ t3_n123 ][4] = t_e4_13;
    data_tetrahedron_edges[ t3_n123 ][5] = t_e5_23;
    
    data_tetrahedron_faces[ t0_n012 ][0] = f0_n01;
    data_tetrahedron_faces[ t0_n012 ][1] = f1_n02;
    data_tetrahedron_faces[ t0_n012 ][2] = f3_n12;
    data_tetrahedron_faces[ t0_n012 ][3] = t_f0_012;
    
    data_tetrahedron_faces[ t1_n013 ][0] = f0_n01;
    data_tetrahedron_faces[ t1_n013 ][1] = f2_n03;
    data_tetrahedron_faces[ t1_n013 ][2] = f4_n13;
    data_tetrahedron_faces[ t1_n013 ][3] = t_f1_013;
    
    data_tetrahedron_faces[ t2_n023 ][0] = f1_n02;
    data_tetrahedron_faces[ t2_n023 ][1] = f2_n03;
    data_tetrahedron_faces[ t2_n023 ][2] = f5_n23;
    data_tetrahedron_faces[ t2_n023 ][3] = t_f2_023;
    
    data_tetrahedron_faces[ t3_n123 ][0] = f3_n12;
    data_tetrahedron_faces[ t3_n123 ][1] = f4_n13;
    data_tetrahedron_faces[ t3_n123 ][2] = f5_n23;
    data_tetrahedron_faces[ t3_n123 ][3] = t_f3_123;
    
    data_face_vertices[ f0_n01 ][0] = vn;
    data_face_vertices[ f0_n01 ][1] = t_v0;
    data_face_vertices[ f0_n01 ][2] = t_v1;
    
    data_face_vertices[ f1_n02 ][0] = vn;
    data_face_vertices[ f1_n02 ][1] = t_v0;
    data_face_vertices[ f1_n02 ][2] = t_v2;
    
    data_face_vertices[ f2_n03 ][0] = vn;
    data_face_vertices[ f2_n03 ][1] = t_v0;
    data_face_vertices[ f2_n03 ][2] = t_v3;
    
    data_face_vertices[ f3_n12 ][0] = vn;
    data_face_vertices[ f3_n12 ][1] = t_v1;
    data_face_vertices[ f3_n12 ][2] = t_v2;
    
    data_face_vertices[ f4_n13 ][0] = vn;
    data_face_vertices[ f4_n13 ][1] = t_v1;
    data_face_vertices[ f4_n13 ][2] = t_v3;
    
    data_face_vertices[ f5_n23 ][0] = vn;
    data_face_vertices[ f5_n23 ][1] = t_v2;
    data_face_vertices[ f5_n23 ][2] = t_v3;
    
    data_face_edges[ f0_n01 ][0] = e0_n0;
    data_face_edges[ f0_n01 ][1] = e1_n1;
    data_face_edges[ f0_n01 ][2] = t_e0_01;
    
    data_face_edges[ f1_n02 ][0] = e0_n0;
    data_face_edges[ f1_n02 ][1] = e2_n2;
    data_face_edges[ f1_n02 ][2] = t_e1_02;
    
    data_face_edges[ f2_n03 ][0] = e0_n0;
    data_face_edges[ f2_n03 ][1] = e3_n3;
    data_face_edges[ f2_n03 ][2] = t_e2_03;
    
    data_face_edges[ f3_n12 ][0] = e1_n1;
    data_face_edges[ f3_n12 ][1] = e2_n2;
    data_face_edges[ f3_n12 ][2] = t_e3_12;
    
    data_face_edges[ f4_n13 ][0] = e1_n1;
    data_face_edges[ f4_n13 ][1] = e3_n3;
    data_face_edges[ f4_n13 ][2] = t_e4_13;
    
    data_face_edges[ f5_n23 ][0] = e2_n2;
    data_face_edges[ f5_n23 ][1] = e3_n3;
    data_face_edges[ f5_n23 ][2] = t_e5_23;
    
    data_edge_vertices[ e0_n0 ][0] = vn;
    data_edge_vertices[ e0_n0 ][1] = t_v0;
    
    data_edge_vertices[ e1_n1 ][0] = vn;
    data_edge_vertices[ e1_n1 ][1] = t_v1;
    
    data_edge_vertices[ e2_n2 ][0] = vn;
    data_edge_vertices[ e2_n2 ][1] = t_v2;
    
    data_edge_vertices[ e3_n3 ][0] = vn;
    data_edge_vertices[ e3_n3 ][1] = t_v3;


    // PROOFREADING DONE UNTIL HERE 23.03.2025
    
    // 4. Check that the first parents are set and actual parents
    //    Higher dimensional parents first, higher dimensional children first
    // 5. Check that the next parents are actual parents
    //    Lower dimensional children first, lower dimensional parents first
    // 6. Check that all parents are listed
    //    Higher dimensional parents first, higher dimensional children first

    /* We fill in the nextparents field for the new tetrahedra, using new indices and recycling old indices */
    /* For the tetrahedron parents of the old faces/edges/vertices, replace t with t0, t1, t2, t3 */
    
    
    data_tetrahedron_nextparents_of_faces[ t0_n012 ][0] = t1_n013;    /* n01 */  data_face_firstparent_tetrahedron[f0_n01] = t0_n012;
    data_tetrahedron_nextparents_of_faces[ t0_n012 ][1] = t2_n023;    /* n02 */  data_face_firstparent_tetrahedron[f1_n02] = t0_n012;
    data_tetrahedron_nextparents_of_faces[ t0_n012 ][2] = t3_n123;    /* n12 */  data_face_firstparent_tetrahedron[f3_n12] = t0_n012;
    data_tetrahedron_nextparents_of_faces[ t0_n012 ][3] = x_t_f0_012; /* 012 */
    
    data_tetrahedron_nextparents_of_faces[ t1_n013 ][0] = nullindex;  /* n01 */
    data_tetrahedron_nextparents_of_faces[ t1_n013 ][1] = t2_n023;    /* n03 */  data_face_firstparent_tetrahedron[f2_n03] = t1_n013;
    data_tetrahedron_nextparents_of_faces[ t1_n013 ][2] = t3_n123;    /* n13 */  data_face_firstparent_tetrahedron[f4_n13] = t1_n013;
    data_tetrahedron_nextparents_of_faces[ t1_n013 ][3] = x_t_f1_013; /* 013 */
    
    data_tetrahedron_nextparents_of_faces[ t2_n023 ][0] = nullindex;  /* n02 */
    data_tetrahedron_nextparents_of_faces[ t2_n023 ][1] = nullindex;  /* n03 */
    data_tetrahedron_nextparents_of_faces[ t2_n023 ][2] = t3_n123;    /* n23 */  data_face_firstparent_tetrahedron[f5_n23] = t2_n023;
    data_tetrahedron_nextparents_of_faces[ t2_n023 ][3] = x_t_f2_023; /* 023 */
    
    data_tetrahedron_nextparents_of_faces[ t3_n123 ][0] = nullindex;  /* n12 */
    data_tetrahedron_nextparents_of_faces[ t3_n123 ][1] = nullindex;  /* n13 */
    data_tetrahedron_nextparents_of_faces[ t3_n123 ][2] = nullindex;  /* n23 */
    data_tetrahedron_nextparents_of_faces[ t3_n123 ][3] = x_t_f3_123; /* 123 */

    { 
      constexpr int N = 4;
      const int faces_to_relink[N] = { t_f0_012, t_f1_013, t_f2_023, t_f3_123 }; 
      const int new_parent_tets[N] = {  t0_n012,  t1_n013,  t2_n023,  t3_n123 }; 
      
      for( int i = 0; i < N; i++ )
      {
        const auto face_to_relink = faces_to_relink[i];
        const auto new_parent_tet = new_parent_tets[i];
        
        if( data_face_firstparent_tetrahedron[ face_to_relink ] == t ) 
          data_face_firstparent_tetrahedron[ face_to_relink ] = new_parent_tet;
        else {
          int current_tetrahedron = data_face_firstparent_tetrahedron[ face_to_relink ];
          while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ] != t 
                &&
                data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ] != nullindex )
            current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ];
          assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ] != nullindex );
          assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ] == t );
          data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, face_to_relink ) ] = new_parent_tet;
        }

        assert( data_face_firstparent_tetrahedron[ face_to_relink ] != nullindex );
      }
    }
    
    
    
    
    
    
    
    
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][0] = t1_n013;   /*   e0_n0 */ data_edge_firstparent_tetrahedron[e0_n0] = t0_n012;
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][1] = t1_n013;   /*   e1_n1 */ data_edge_firstparent_tetrahedron[e1_n1] = t0_n012;
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][2] = t2_n023;   /*   e2_n2 */ data_edge_firstparent_tetrahedron[e2_n2] = t0_n012;
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][3] = t1_n013;   /* t_e0_01 */
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][4] = t2_n023;   /* t_e1_02 */ 
    data_tetrahedron_nextparents_of_edges[ t0_n012 ][5] = t3_n123;   /* t_e3_12 */
    
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][0] = t2_n023;   /*   e0_n0 */
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][1] = t3_n123;   /*   e1_n1 */
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][2] = t2_n023;   /*   e3_n3 */ data_edge_firstparent_tetrahedron[e3_n3] = t1_n013;
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][3] = x_t_e0_01; /* t_e0_01 */
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][4] = t2_n023;   /* t_e2_03 */
    data_tetrahedron_nextparents_of_edges[ t1_n013 ][5] = t3_n123;   /* t_e4_13 */
    
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][0] = nullindex; /*   e0_n0 */
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][1] = t3_n123;   /*   e2_n2 */
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][2] = t3_n123;   /*   e3_n3 */
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][3] = x_t_e1_02; /* t_e1_02 */
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][4] = x_t_e2_03; /* t_e2_03 */
    data_tetrahedron_nextparents_of_edges[ t2_n023 ][5] = t3_n123;   /* t_e5_23 */
    
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][0] = nullindex; /*   e1_n1 */
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][1] = nullindex; /*   e2_n2 */
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][2] = nullindex; /*   e3_n3 */
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][3] = x_t_e3_12; /* t_e3_12 */
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][4] = x_t_e4_13; /* t_e4_13 */
    data_tetrahedron_nextparents_of_edges[ t3_n123 ][5] = x_t_e5_23; /* t_e5_23 */
    
    { 
      constexpr int N = 6;
      const int edges_to_relink[N] = { t_e0_01, t_e1_02, t_e2_03, t_e3_12, t_e4_13, t_e5_23 }; 
      const int new_parent_tets[N] = { t0_n012, t0_n012, t1_n013, t0_n012, t1_n013, t2_n023 }; 
      
      for( int i = 0; i < N; i++ )
      {
        const auto edge_to_relink = edges_to_relink[i];
        const auto new_parent_tet = new_parent_tets[i];
        
        if( data_edge_firstparent_tetrahedron[ edge_to_relink ] == t ) 
          data_edge_firstparent_tetrahedron[ edge_to_relink ] = new_parent_tet;
        else {
          int current_tetrahedron = data_edge_firstparent_tetrahedron[ edge_to_relink ];
          while( data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ] != t 
                &&
                data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ] != nullindex )
            current_tetrahedron = data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ];
          assert( data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ] != nullindex );
          assert( data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ] == t );
          data_tetrahedron_nextparents_of_edges[ current_tetrahedron ][ indexof_tetrahedron_edge( current_tetrahedron, edge_to_relink ) ] = new_parent_tet;
        }

        assert( data_edge_firstparent_tetrahedron[ edge_to_relink ] != nullindex );
      }
    }
    
    
    
    
    
    

    data_tetrahedron_nextparents_of_vertices[ t0_n012 ][0] = t1_n013;   /* n */ data_vertex_firstparent_tetrahedron[vn] = t0_n012;
    data_tetrahedron_nextparents_of_vertices[ t0_n012 ][1] = t1_n013;   /* 0 */
    data_tetrahedron_nextparents_of_vertices[ t0_n012 ][2] = t1_n013;   /* 1 */
    data_tetrahedron_nextparents_of_vertices[ t0_n012 ][3] = t2_n023;   /* 2 */
    
    data_tetrahedron_nextparents_of_vertices[ t1_n013 ][0] = t2_n023;   /* n */
    data_tetrahedron_nextparents_of_vertices[ t1_n013 ][1] = t2_n023;   /* 0 */
    data_tetrahedron_nextparents_of_vertices[ t1_n013 ][2] = t3_n123;   /* 1 */
    data_tetrahedron_nextparents_of_vertices[ t1_n013 ][3] = t2_n023;   /* 3 */
    
    data_tetrahedron_nextparents_of_vertices[ t2_n023 ][0] = t3_n123;   /* n */ 
    data_tetrahedron_nextparents_of_vertices[ t2_n023 ][1] = x_t_v0;    /* 0 */ // data_vertex_firstparent_tetrahedron[ t_v0 ];
    data_tetrahedron_nextparents_of_vertices[ t2_n023 ][2] = t3_n123;   /* 2 */
    data_tetrahedron_nextparents_of_vertices[ t2_n023 ][3] = t3_n123;   /* 3 */
    
    data_tetrahedron_nextparents_of_vertices[ t3_n123 ][0] = nullindex; /* n */
    data_tetrahedron_nextparents_of_vertices[ t3_n123 ][1] = x_t_v1;    /* 1 */ // data_vertex_firstparent_tetrahedron[ t_v1 ];
    data_tetrahedron_nextparents_of_vertices[ t3_n123 ][2] = x_t_v2;    /* 2 */ // data_vertex_firstparent_tetrahedron[ t_v2 ];
    data_tetrahedron_nextparents_of_vertices[ t3_n123 ][3] = x_t_v3;    /* 3 */ // data_vertex_firstparent_tetrahedron[ t_v3 ];
    
    { 
      constexpr int N = 4;
      const int vertices_to_relink[N] = {    t_v0,    t_v1,    t_v2,    t_v3 }; 
      const int    new_parent_tets[N] = { t0_n012, t0_n012, t0_n012, t1_n013 }; 
      
      for( int i = 0; i < N; i++ )
      {
        const auto vertex_to_relink = vertices_to_relink[i];
        const auto new_parent_tet   =    new_parent_tets[i];
        
        if( data_vertex_firstparent_tetrahedron[ vertex_to_relink ] == t ) 
          data_vertex_firstparent_tetrahedron[ vertex_to_relink ] = new_parent_tet;
        else {
          int current_tetrahedron = data_vertex_firstparent_tetrahedron[ vertex_to_relink ];
          assert( current_tetrahedron != t );
          assert( data_tetrahedron_vertices[current_tetrahedron][0] == vertex_to_relink || data_tetrahedron_vertices[current_tetrahedron][1] == vertex_to_relink || data_tetrahedron_vertices[current_tetrahedron][2] == vertex_to_relink || data_tetrahedron_vertices[current_tetrahedron][3] == vertex_to_relink ); 
          while( data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ] != t 
                &&
                data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ] != nullindex )
            current_tetrahedron = data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ];
          assert( data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ] != nullindex );
          assert( data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ] == t );
          data_tetrahedron_nextparents_of_vertices[ current_tetrahedron ][ indexof_tetrahedron_vertex( current_tetrahedron, vertex_to_relink ) ] = new_parent_tet;
        }

        assert( data_vertex_firstparent_tetrahedron[ vertex_to_relink ] != nullindex );
      }
    }
    
    
    
    /* TODO: we fill in the nextparents field for the new subsimplices, using only new indices and finish with nullindex         */
    /*       with regard to old faces/edges/vertices, we change the first parent and append the old list, as a post-modification */
    
    data_face_nextparents_of_edges[ f0_n01 ][0] = f1_n02;     /*   e0_n0 */ data_edge_firstparent_face[   e0_n0 ] = f0_n01;
    data_face_nextparents_of_edges[ f0_n01 ][1] = f3_n12;     /*   e1_n1 */ data_edge_firstparent_face[   e1_n1 ] = f0_n01;
    data_face_nextparents_of_edges[ f0_n01 ][2] = fp_e0_01_f; /* t_e0_01 */ data_edge_firstparent_face[ t_e0_01 ] = f0_n01;
    
    data_face_nextparents_of_edges[ f1_n02 ][0] = f2_n03;     /*   e0_n0 */
    data_face_nextparents_of_edges[ f1_n02 ][1] = f3_n12;     /*   e2_n2 */ data_edge_firstparent_face[   e2_n2 ] = f1_n02;
    data_face_nextparents_of_edges[ f1_n02 ][2] = fp_e1_02_f; /* t_e1_02 */ data_edge_firstparent_face[ t_e1_02 ] = f1_n02;
    
    data_face_nextparents_of_edges[ f2_n03 ][0] = nullindex;  /*   e0_n0 */
    data_face_nextparents_of_edges[ f2_n03 ][1] = f4_n13;     /*   e3_n3 */ data_edge_firstparent_face[   e3_n3 ] = f2_n03;
    data_face_nextparents_of_edges[ f2_n03 ][2] = fp_e2_03_f; /* t_e2_03 */ data_edge_firstparent_face[ t_e2_03 ] = f2_n03;
    
    data_face_nextparents_of_edges[ f3_n12 ][0] = f4_n13;     /*   e1_n1 */
    data_face_nextparents_of_edges[ f3_n12 ][1] = f5_n23;     /*   e2_n2 */
    data_face_nextparents_of_edges[ f3_n12 ][2] = fp_e3_12_f; /* t_e3_12 */ data_edge_firstparent_face[ t_e3_12 ] = f3_n12;
    
    data_face_nextparents_of_edges[ f4_n13 ][0] = nullindex;  /*   e1_n1 */
    data_face_nextparents_of_edges[ f4_n13 ][1] = f5_n23;     /*   e3_n3 */
    data_face_nextparents_of_edges[ f4_n13 ][2] = fp_e4_13_f; /* t_e4_13 */ data_edge_firstparent_face[ t_e4_13 ] = f4_n13;
    
    data_face_nextparents_of_edges[ f5_n23 ][0] = nullindex;  /*   e2_n2 */
    data_face_nextparents_of_edges[ f5_n23 ][1] = nullindex;  /*   e3_n3 */
    data_face_nextparents_of_edges[ f5_n23 ][2] = fp_e5_23_f; /* t_e5_23 */ data_edge_firstparent_face[ t_e5_23 ] = f5_n23;
    
    
    data_face_nextparents_of_vertices[ f0_n01 ][0] = f1_n02;    /* n */  data_vertex_firstparent_face[  vn] = f0_n01;
    data_face_nextparents_of_vertices[ f0_n01 ][1] = f1_n02;    /* 0 */  data_vertex_firstparent_face[t_v0] = f0_n01;
    data_face_nextparents_of_vertices[ f0_n01 ][2] = f3_n12;    /* 1 */  data_vertex_firstparent_face[t_v1] = f0_n01;
    
    data_face_nextparents_of_vertices[ f1_n02 ][0] = f2_n03;    /* n */
    data_face_nextparents_of_vertices[ f1_n02 ][1] = f2_n03;    /* 0 */
    data_face_nextparents_of_vertices[ f1_n02 ][2] = f3_n12;    /* 2 */  data_vertex_firstparent_face[t_v2] = f1_n02;
    
    data_face_nextparents_of_vertices[ f2_n03 ][0] = f3_n12;    /* n */
    data_face_nextparents_of_vertices[ f2_n03 ][1] = fp_v0_f;   /* 0 */  
    data_face_nextparents_of_vertices[ f2_n03 ][2] = f4_n13;    /* 3 */  data_vertex_firstparent_face[t_v3] = f2_n03;
    
    data_face_nextparents_of_vertices[ f3_n12 ][0] = f4_n13;    /* n */
    data_face_nextparents_of_vertices[ f3_n12 ][1] = f4_n13;    /* 1 */
    data_face_nextparents_of_vertices[ f3_n12 ][2] = f5_n23;    /* 2 */
    
    data_face_nextparents_of_vertices[ f4_n13 ][0] = f5_n23;    /* n */
    data_face_nextparents_of_vertices[ f4_n13 ][1] = fp_v1_f;   /* 1 */
    data_face_nextparents_of_vertices[ f4_n13 ][2] = f5_n23;    /* 3 */
    
    data_face_nextparents_of_vertices[ f5_n23 ][0] = nullindex; /* n */
    data_face_nextparents_of_vertices[ f5_n23 ][1] = fp_v2_f;   /* 2 */
    data_face_nextparents_of_vertices[ f5_n23 ][2] = fp_v3_f;   /* 3 */
    
    
    data_edge_nextparents_of_vertices[ e0_n0 ][0] = e1_n1;     /* n */  data_vertex_firstparent_edge[  vn] = e0_n0;
    data_edge_nextparents_of_vertices[ e0_n0 ][1] = fp_v0_e;   /* 0 */  data_vertex_firstparent_edge[t_v0] = e0_n0;
    
    data_edge_nextparents_of_vertices[ e1_n1 ][0] = e2_n2;     /* n */
    data_edge_nextparents_of_vertices[ e1_n1 ][1] = fp_v1_e;   /* 1 */  data_vertex_firstparent_edge[t_v1] = e1_n1;
    
    data_edge_nextparents_of_vertices[ e2_n2 ][0] = e3_n3;     /* n */
    data_edge_nextparents_of_vertices[ e2_n2 ][1] = fp_v2_e;   /* 2 */  data_vertex_firstparent_edge[t_v2] = e2_n2;
    
    data_edge_nextparents_of_vertices[ e3_n3 ][0] = nullindex; /* n */
    data_edge_nextparents_of_vertices[ e3_n3 ][1] = fp_v3_e;   /* 3 */  data_vertex_firstparent_edge[t_v3] = e3_n3;
    
    
    
    /* Counters */
    counter_tetrahedra = counter_tetrahedra + 3;
    counter_faces      = counter_faces      + 6;
    counter_edges      = counter_edges      + 4;
    counter_vertices   = counter_vertices   + 1;
    
    /* DONE */
    
    check();
}


void MeshSimplicial3D::midpoint_refinement_global()
{
    check();
    
    int N = counter_tetrahedra;
    
    // adjust capacity of containers for efficiency reasons

    data_tetrahedron_faces.reserve               ( counter_tetrahedra + 3 * N );
    data_face_firstparent_tetrahedron.reserve    ( counter_faces      + 6 * N );
    data_tetrahedron_nextparents_of_faces.reserve( counter_tetrahedra + 3 * N );
    
    data_tetrahedron_edges.reserve               ( counter_tetrahedra + 3 * N );
    data_edge_firstparent_tetrahedron.reserve    ( counter_edges      + 4 * N );
    data_tetrahedron_nextparents_of_edges.reserve( counter_tetrahedra + 3 * N );
    
    data_tetrahedron_vertices.reserve               ( counter_tetrahedra + 3 * N );
    data_vertex_firstparent_tetrahedron.reserve     ( counter_vertices   + 1 * N );
    data_tetrahedron_nextparents_of_vertices.reserve( counter_tetrahedra + 3 * N );
    
    data_face_edges.reserve               ( counter_faces + 6 * N );
    data_edge_firstparent_face.reserve    ( counter_edges + 4 * N );
    data_face_nextparents_of_edges.reserve( counter_faces + 6 * N );
    
    data_face_vertices.reserve               ( counter_faces    + 6 * N );
    data_vertex_firstparent_face.reserve     ( counter_vertices + 1 * N );
    data_face_nextparents_of_vertices.reserve( counter_faces    + 6 * N );
    
    data_edge_vertices.reserve               ( counter_edges    + 4 * N );
    data_vertex_firstparent_edge.reserve     ( counter_vertices + 1 * N );
    data_edge_nextparents_of_vertices.reserve( counter_edges    + 4 * N );

    flags_tetrahedra.reserve( counter_tetrahedra + 3 * N );
    flags_faces.reserve     ( counter_faces      + 6 * N );
    flags_edges.reserve     ( counter_edges      + 4 * N );
    flags_vertices.reserve  ( counter_vertices   + 1 * N );



    for( int t = 0; t < N; t++ ) {
      
      midpoint_refinement( t );
      
    }
    
    check();
}











