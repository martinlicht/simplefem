
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
#include "simplicialmesh.hpp"
#include "io.coordinates.hpp"




void writeCoordinates( const char* filename, const Coordinates& coords, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeCoordinates( myfile, coords, sugar );
    myfile.close();
}

Coordinates readCoordinatesPath( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    Coordinates coords = readCoordinates( myfile );
    myfile.close();
    return coords;
}


void writeCoordinates( std::ostream& out, const Coordinates& coords, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing coordinates" << std::endl;
    if( sugar ) out << "Dimension of coordinates: " << std::endl;
    out << coords.getdimension();
    if( sugar ) out << "Number of points: " << std::endl;
    out << coords.getnumber();
    
    /* data */
    for( int p = 0; p < coords.getnumber(); p++ ) {
        if( sugar ) out << p << ": ";
        for( int d = 0; d < coords.getdimension(); d++ )
            out << coords.getdata(p,d) << " ";
        out << std::endl;
    }
}

Coordinates readCoordinatesPath( std::istream& in )
{
    int dimension;
    int numpoints;
    in >> dimension >> numpoints;
    Coordinates coords( dimension, numpoints );
    for( int p = 0; p < numpoints; p++ ) {
        Float tempvalue;
        for( int d = 0; d < dimension; d++ ) {
            Float tempvalue;
            in >> tempvalue;
            coords.setdata( p, d, tempvalue );
        }
    }
    return coords;
}



        
        
        
        
        