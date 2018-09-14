#ifndef INCLUDEGUARD_EXAMPLES_2D
#define INCLUDEGUARD_EXAMPLES_2D


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "mesh.hpp"
#include "mesh.simplicial2D.hpp"







inline MeshSimplicial2D UnitTriangle2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 3, {
        0., 0., // 0
        1., 0., // 1
        0., 1.  // 2
      } ),
      {
        { 0, 1, 2 }
      }
    );
}

inline MeshSimplicial2D UnitSquare2D()
{
//     return MeshSimplicial2D(
//       2,
//       Coordinates( 2, 4, {
//         -1., -1., // 0
//         -1.,  1., // 1
//          1., -1., // 2
//          1.,  1.  // 3
//       } ),
//       {
//         { 0, 1, 3 },
//         { 0, 2, 3 }
//       }
//     );
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 4, {
        -1.,  1., // 0
         1., -1., // 1
        -1., -1., // 2
         1.,  1.  // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 }
      }
    );
}

inline MeshSimplicial2D StandardSquare2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 9, {
         0.,  0., // 0
         1.,  0., // 1
         1.,  1., // 2
         0.,  1., // 3
        -1.,  1., // 4
        -1.,  0., // 5
        -1., -1., // 6
         0., -1., // 7
         1., -1.  // 8
      } ),
      {
        { 0, 1, 2 },
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 6 },
        { 0, 6, 7 },
        { 0, 7, 8 },
        { 0, 8, 1 }
      }
    );
}



inline MeshSimplicial2D LShapedDomain2D()
{
    return MeshSimplicial2D(
      3,
      Coordinates( 3, 8, {
         0.,  0.,  0., // 0
         1.,  0.,  0., // 1
         1.,  1.,  0., // 2
         0.,  1.,  0., // 3
        -1.,  1.,  0., // 4
        -1.,  0.,  0., // 5
        -1., -1.,  0., // 6
         0., -1.,  0.  // 7
      } ),
      {
        { 0, 1, 2 },
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 6 },
        { 0, 6, 7 }
      }
    );
}



inline MeshSimplicial2D SlitDomain2D()
{
    return MeshSimplicial2D(
      3,
      Coordinates( 3, 10, {
         0.,  0., 0, // 0
         1.,  0., 0.5, // 1
         1.,  1., 0, // 2
         0.,  1., 0, // 3
        -1.,  1., 0, // 4
        -1.,  0., 0, // 5
        -1., -1., 0, // 6
         0., -1., 0, // 7
         1., -1., 0, // 8
         1.,  0., -0.5,  // 9
      } ),
      {
        { 0, 1, 2 }, // 0
        { 0, 2, 3 }, // 1
        { 0, 3, 4 }, // 2
        { 0, 4, 5 }, // 3
        { 0, 5, 6 }, // 4
        { 0, 6, 7 }, // 5
        { 0, 7, 8 }, // 6
        { 0, 8, 9 }  // 7
      }
    );
}


inline MeshSimplicial2D TetrahedralSurface2D()
{
    return MeshSimplicial2D(
      3,
      Coordinates( 3, 4, {
         0.,  0.,  0., // 0
         1.,  0.,  0., // 1
         0.,  1.,  0., // 2
         0.,  0.,  1.,  // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
      }
    );
}










#endif
