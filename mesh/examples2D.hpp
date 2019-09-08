#ifndef INCLUDEGUARD_EXAMPLES_2D
#define INCLUDEGUARD_EXAMPLES_2D


#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>


#include "mesh.hpp"
#include "mesh.simplicial2D.hpp"







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




inline MeshSimplicial2D UnitSquare2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 4, {
        -1., -1., // 0
        -1.,  1., // 1
         1., -1., // 2
         1.,  1.  // 3
      } ),
      {
        { 0, 1, 3 },
        { 0, 2, 3 }
      }
    );
}

inline MeshSimplicial2D UnitSquare2D_Alternative()
{
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










inline MeshSimplicial2D UnitDisk( int L = 1 )
{
    assert( L >= 1 );
    
    // 1. Calculate the number of vertices and triangles 
    
    //  3 + 6 + 12 + 24 + ... 
    //= 3 ( 1 + 2 + 4 + 8 + .... )
    //= 3 ( 2^( L+1 ) - 1 )
    int num_vertices = 3 * ((1<<L)-1);
    
    //  1 + (3+6) + (6+12) + (12+24) + ...
    //= 1 +   9   +    18  +     36  + ...
    int num_triangles = 9 * (1<<(L-1)) - 8;
    
    
    
    // 2. Create the coordinates 
    
    std::vector<Float> coords;
    coords.reserve( 2 * num_vertices );
    
    // 2.1 Calculate the radii
    
    std::vector<Float> Rs(L);
    Rs[0] = 1.;
    for( int l = 1; l < L; l++ ) Rs[l] = Rs[l-1] + 1.5 * 2 * 3.14159 * Rs[l-1] / ( 3. * (1<<(l)) );
    Float Rmax = *std::max_element(Rs.begin(),Rs.end());
    for( Float& R : Rs ) R /= Rmax;
    
    // 2.2 fill in the values
    
    for( int l = 1; l <= L; l++ ) {
        
        int N = 3 * (1<<(l-1));
        
        Float radius = Rs[l-1];
        
        for( int a = 0; a < N; a++ )
        {
            coords.push_back( radius * std::cos( 2 * 3.14159 * a / (Float)N ) );
            coords.push_back( radius * std::sin( 2 * 3.14159 * a / (Float)N ) );
        }
    }
    
    std::cout << coords.size() / 2 << space << num_vertices << nl;
    assert( coords.size() == 2 * num_vertices );
    
    std::vector<std::array<int,3>> tris;
    tris.reserve( num_triangles );
    
    tris.push_back( {0,1,2} );
    
    for( int l = 1; l < L; l++ )
    {
        int base_inner  = 3 * ( integerpower( 2, l-1 ) - 1 );
        int count_inner = 3 * integerpower( 2, l-1 );
        
        int base_outer  = base_inner + count_inner;
        int count_outer = 2 * count_inner;
        
        for( int i = 0; i < count_inner; i++ )
        {
            
            tris.push_back( { 
                base_inner + i, 
                base_outer + 2*i,
                base_outer + (2*i + 1)%count_outer,
                } );
            tris.push_back( { 
                base_inner + i, 
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            tris.push_back( { 
                base_outer + (2*i+2)%count_outer,
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            
        }
    }
    
    std::cout << coords.size() / 2 << space << num_vertices << nl;
    std::cout << tris.size() << space << num_triangles << nl;
    assert( tris.size() == num_triangles );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         std::cout << v << space; std::cout << nl; }
    
    for( auto& t : tris ) std::sort( t.begin(), t.end() );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         std::cout << v << space; std::cout << nl; }
    
    return MeshSimplicial2D(
      2,
      Coordinates( 2, num_vertices, coords ),
      tris
    );
}









inline MeshSimplicial2D Annulus( int Linner, int Louter = 1 )
{
    assert( Linner >= 1 && Louter > Linner );
    
    // 1. Calculate the number of vertices and triangles 
    
    //  3 + 6 + 12 + 24 + ... 
    //= 3 ( 1 + 2 + 4 + 8 + .... )
    //= 3 ( 2^( L+1 ) - 1 )
    int num_vertices  = 3 * ((1<<Louter)-1)     - ( 3 * ((1<<(Linner-1))-1) );
    
    //  1 + (3+6) + (6+12) + (12+24) + ...
    //= 1 +   9   +    18  +     36  + ...
    int num_triangles = 9 * (1<<(Louter-1)) - 8 - ( 9 * (1<<(Linner-1)) - 8 );
    
    
    
    // 2. Create the coordinates 
    
    std::vector<Float> coords;
    coords.reserve( 2 * num_vertices );
    
    // 2.1 Calculate the radii
    
    std::vector<Float> Rs(Louter);
    for( int l = 0; l < Linner; l++ ) Rs[l] = 1.;
    for( int l = Linner; l < Louter; l++ ) Rs[l] = Rs[l-1] + 1.5 * 2 * 3.14159 * Rs[l-1] / ( 3. * (1<<(l)) );
    Float Rmax = *std::max_element(Rs.begin(),Rs.end());
    for( Float& R : Rs ) R = pow( R / Rmax, 1.5 );
    
    // 2.2 fill in the values
    
    for( int l = Linner; l <= Louter; l++ ) {
        
        int N = 3 * (1<<(l-1));
        
        Float radius = Rs[l-1];
        
        for( int a = 0; a < N; a++ )
        {
            coords.push_back( radius * std::cos( 2 * 3.14159 * a / (Float)N ) );
            coords.push_back( radius * std::sin( 2 * 3.14159 * a / (Float)N ) );
        }
    }
    
    std::cout << coords.size() / 2 << space << num_vertices << nl;
    assert( coords.size() == 2 * num_vertices );
    
    std::vector<std::array<int,3>> tris;
    tris.reserve( num_triangles );
    
    for( int l = Linner; l < Louter; l++ )
    {
        int base_inner  = 3 * ( integerpower( 2, l-1 ) - 1 ) - 3 * ( integerpower( 2, Linner-1 ) - 1 );
        int count_inner = 3 * integerpower( 2, l-1 );
        
        int base_outer  = base_inner + count_inner;
        int count_outer = 2 * count_inner;
        
        for( int i = 0; i < count_inner; i++ )
        {
            
            tris.push_back( { 
                base_inner + i, 
                base_outer + 2*i,
                base_outer + (2*i + 1)%count_outer,
                } );
            tris.push_back( { 
                base_inner + i, 
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            tris.push_back( { 
                base_outer + (2*i+2)%count_outer,
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            
        }
    }
    
    std::cout << coords.size() / 2 << space << num_vertices << nl;
    std::cout << tris.size() << space << num_triangles << nl;
    assert( tris.size() == num_triangles );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         std::cout << v << space; std::cout << nl; }
    
    for( auto& t : tris ) std::sort( t.begin(), t.end() );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         std::cout << v << space; std::cout << nl; }
    
    return MeshSimplicial2D(
      2,
      Coordinates( 2, num_vertices, coords ),
      tris
    );
}





inline MeshSimplicial2D SphericalSurface2D( int L = 0 )
{
    
    MeshSimplicial2D ret(
      3,
      Coordinates( 3, 4, {
          sqrt(8./9.),          0.0, -1./3., // 0
         -sqrt(2./9.),  sqrt(2./3.), -1./3., // 1
         -sqrt(2./9.), -sqrt(2./3.), -1./3., // 2
                  0.0,          0.0,    1.0, // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
      }
    );
    
    for( int l = 0; l < L; l++ )
        ret.uniformrefinement();
    
    for( int n = 0; n < ret.getcoordinates().getnumber(); n++ )
    {
        FloatVector point = ret.getcoordinates().getvectorclone( n );
        point.normalize();
        ret.getcoordinates().loadvector( n, point );
    }
    
    return ret;
    
}






#endif
