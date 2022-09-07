#ifndef INCLUDEGUARD_EXAMPLES_3D_HPP
#define INCLUDEGUARD_EXAMPLES_3D_HPP


#include "mesh.hpp"
#include "mesh.simplicial3D.hpp"






inline MeshSimplicial3D UnitCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         0., 0., 0., // 000
         0., 0., 1., // 001
         0., 1., 0., // 010
         0., 1., 1., // 011
         1., 0., 0., // 100
         1., 0., 1., // 101
         1., 1., 0., // 110
         1., 1., 1.  // 111
      } ),
      {
        { 0b000, 0b100, 0b110, 0b111 }, 
        { 0b000, 0b100, 0b101, 0b111 },
        { 0b000, 0b010, 0b110, 0b111 },
        { 0b000, 0b010, 0b011, 0b111 },
        { 0b000, 0b001, 0b101, 0b111 },
        { 0b000, 0b001, 0b011, 0b111 }
      }
    );
}

inline MeshSimplicial3D StandardCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         -1., -1., -1., // 000
         -1., -1., +1., // 001
         -1., +1., -1., // 010
         -1., +1., +1., // 011
         +1., -1., -1., // 100
         +1., -1., +1., // 101
         +1., +1., -1., // 110
         +1., +1., +1.  // 111
      } ),
      {
        { 0b000, 0b100, 0b110, 0b111 }, 
        { 0b000, 0b100, 0b101, 0b111 },
        { 0b000, 0b010, 0b110, 0b111 },
        { 0b000, 0b010, 0b011, 0b111 },
        { 0b000, 0b001, 0b101, 0b111 },
        { 0b000, 0b001, 0b011, 0b111 }
      }
    );
}






inline MeshSimplicial3D OctoDiamond3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         1., 0., 0., // 1
         0., 1., 0., // 2
        -1., 0., 0., // 3
         0.,-1., 0., // 4
         0., 0., 1., // 5
         0., 0.,-1.  // 6
      } ),
      {
        { 0, 1, 2, 5 },
        { 0, 2, 3, 5 },
        { 0, 3, 4, 5 },
        { 0, 1, 4, 5 },
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 },
        { 0, 3, 4, 6 },
        { 0, 1, 4, 6 } 
      }
    );
}




inline MeshSimplicial3D UnitSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         0., 0., 1., // 1
         0., 1., 0., // 2
         1., 0., 0.  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}

inline MeshSimplicial3D SomeSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0.2, 0.0,-0.1, // 0
         0.0, 0.3, 0.9, // 1
         0.1, 1.2, 0.1, // 2
         1.5,-0.1,-0.2  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}

inline MeshSimplicial3D RegularSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         1., 1., 0., // 1
         1., 0., 1., // 2
         0., 1., 1.  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}










#endif
