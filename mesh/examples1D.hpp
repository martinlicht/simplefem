#ifndef INCLUDEGUARD_EXAMPLES_1D
#define INCLUDEGUARD_EXAMPLES_1D


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "mesh.hpp"
#include "mesh.simplicial1D.hpp"







inline MeshSimplicial1D UnitSquare1D()
{
    return MeshSimplicial1D(
      2,
      Coordinates( 2, 2, {
        -1., 0.333, // 0
         1., 0.444  // 1
      } ),
      {
        { 0, 1 }
      }
    );
}

inline MeshSimplicial1D StandardSquare1D()
{
    return MeshSimplicial1D(
      1,
      Coordinates( 1, 2, {
         0., // 0
         1., // 1
      } ),
      {
        { 0, 1 }
      }
    );
}











#endif
