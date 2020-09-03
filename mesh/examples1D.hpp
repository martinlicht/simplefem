#ifndef INCLUDEGUARD_EXAMPLES_1D_HPP
#define INCLUDEGUARD_EXAMPLES_1D_HPP


#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>


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
