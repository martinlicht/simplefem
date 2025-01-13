#ifndef INCLUDEGUARD_EXAMPLES_ND_HPP
#define INCLUDEGUARD_EXAMPLES_ND_HPP


#include "mesh.hpp"
#include "mesh.simplicialND.hpp"




inline MeshSimplicialND UnitSquareND()
{
    return MeshSimplicialND(
      1, 3, 
      Coordinates( 3, 2, {
        -1., 0., 1., // 0
         1., 3., 2., // 1
      } ),
      {
        0, 1 
      }
    );
}



inline MeshSimplicialND StandardSquare3D()
{
    return MeshSimplicialND(
      3, 3, 
      Coordinates( 3, 8, {
         0., 0., 0., // 0
         0., 0., 1., // 0
         0., 1., 0., // 0
         0., 1., 1., // 0
         1., 0., 0., // 0
         1., 0., 1., // 0
         1., 1., 0., // 0
         1., 1., 1.  // 0
      } ),
      {
        0b000, 0b100, 0b110, 0b111, 
        0b000, 0b100, 0b101, 0b111,
        0b000, 0b010, 0b110, 0b111,
        0b000, 0b010, 0b011, 0b111,
        0b000, 0b001, 0b101, 0b111,
        0b000, 0b001, 0b011, 0b111 
      }
    );
}



inline MeshSimplicialND HypertetrahedralSurface4D()
{
    return MeshSimplicialND(
      3, 4,
      Coordinates( 4, 5, {
         0.,  0.,  0., 0., // 0
         1.,  0.,  0., 0., // 1
         0.,  1.,  0., 0., // 2
         0.,  0.,  1., 0., // 3
         0.,  0.,  0., 1., // 4
      } ),
      {
        0, 1, 2, 3,
        0, 1, 2, 4,
        0, 1, 3, 4,
        0, 2, 3, 4,
        1, 2, 3, 4 
      }
    );
}








#endif
