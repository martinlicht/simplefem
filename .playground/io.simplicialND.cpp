
#include <algorithm>
#include <fstream>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <utility>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "mesh.simplicialND.hpp"
#include "io.simplicialND.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicialND( const char* filename, const MeshSimplicialND& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicialND( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicialND readMeshSimplicialND( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicialND mesh = readMeshSimplicialND( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicialND( std::ostream& out, const MeshSimplicialND& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial ND Mesh..." << nl;
    
    if( sugar ) out << "inner dimension: " << nl;;
    out << mesh.getinnerdimension() << nl;
    
    if( sugar ) out << "number of top-dimensional simplices: " << nl;;
    out << mesh.count_simplices( mesh.getinnerdimension() ) << nl;
    
    if( sugar ) out << "number of vertices: " << nl;
    out << mesh.count_simplices(0) << nl;
    
    if( sugar ) out << "external dimension: " << nl;
    out << mesh.getcoordinates().getdimension() << nl;
    
    /* simplices -> vertices */
    if( sugar ) out << "for each simplex, the vertices: " << nl;
    for( int S = 0; S < mesh.count_simplices( mesh.getinnerdimension() ); S++ ) {
        
        if( sugar ) out << S << ": ";
        
        for( int d = 0; d <= mesh.getinnerdimension(); d++ )
          out << mesh.get_subsimplex( mesh.getinnerdimension(), 0, S, d ) << space;
        
        out << nl;
    }
    
    assert( out.good() );
    
    writeCoordinates( out, mesh.getcoordinates(), sugar );
    
    assert( out.good() );
    
}




MeshSimplicialND readMeshSimplicialND( std::istream& in )
{
    int innerdim, counter_simplices, counter_vertices, outerdim;
    
    in >> innerdim
       >> counter_simplices
       >> counter_vertices
       >> outerdim;
    
    int nullindex = MeshSimplicialND::nullindex;
    
    /* edge -> vertices */
    
    std::vector<int> simplex_vertices( counter_simplices * ( innerdim + 1 ), nullindex );
    
    for( int S = 0; S < counter_simplices; S++ )
      for( int d = 0; d <= innerdim; d++ )
        in >> simplex_vertices[ S * (innerdim+1) +  d ];
    
    assert( in.good() );
    
    /* coordinates */
    
    Coordinates coords = readCoordinates( in );
    
    assert( in.good() );
    
    /* return */
    
    return MeshSimplicialND( innerdim, outerdim, coords, simplex_vertices );
    
}



        
        
        
        
        
