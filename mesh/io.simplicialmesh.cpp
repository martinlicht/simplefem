
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>



#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "coordinates.hpp"
#include "simplicialmesh.hpp"
#include "io.simplicialmesh.hpp"
#include "io.coordinates.hpp"




void writeSimplicialMesh( const char* filename, const SimplicialMesh& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeSimplicialMesh( myfile, mesh, sugar );
    myfile.close();
}

SimplicialMesh readSimplicialMeshPath( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    SimplicialMesh mesh = readSimplicialMesh( myfile );
    myfile.close();
    return mesh;
}


void writeSimplicialMesh( std::ostream& out, const SimplicialMesh& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Mesh output" << std::endl;
    if( sugar ) out << "Inner Dimension: "; 
    out << mesh.getinnerdimension() << std::endl;
    if( sugar ) out << "Outer Dimension: "; 
    out << mesh.getouterdimension() << std::endl;
    
    /* Coordinates */
    writeCoordinates( out, mesh.getcoordinates() );
    
    /* Combinatorial data - count of subsimplices */
    for( int d = 0; d <= mesh.getinnerdimension(); d++ ) {
        out << mesh.countsimplices(d) << " ";
    }
    out << std::endl;
    
    /* Combinatorial data - subsimplices */
    if( sugar ) out << "Sub simplex lists" << std::endl;
    
    const std::map< std::pair<int,int>, std::vector<IndexMap> >& sub = mesh.getsub();
    
    if( sugar ) out << "Number of sub simplex lists: "; 
    out << sub.size() << std::endl;
    
    for( const auto& subentry : sub ) {
    
        if( sugar ) out << "sub simplex list from/to: ";
        out << subentry.first.first << " " << subentry.first.second << std::endl;
        const auto& vecoim = subentry.second;
        if( sugar ) out << "number of lists: ";
        out << vecoim.size();
        for( int supsim = 0; supsim < vecoim.size(); supsim++ ) {
            if( sugar ) out << supsim << ": ";
            int CSS = countsubsimplices( subentry.first.first, subentry.first.second );
            for( int sub = 0; sub < CSS; sub++ ) {
                out << vecoim.at(supsim).at( sub );
            }
            out << std::endl;
        }
    }
    
    /* Combinatorial data - supersimplices */
    if( sugar ) out << "Super simplex lists" << std::endl;
    
    const std::map< std::pair<int,int>, std::vector<std::list<int>> >& super = mesh.getsuper();
    
    if( sugar ) out << "Number of super simplex lists: "; 
    out << super.size() << std::endl;
    
    for( const auto& superentry : super ) {
    
        if( sugar ) out << "super simplex list from/to: ";
        out << superentry.first.first << " " << superentry.first.second << std::endl;
        const auto& vecolist = superentry.second;
        if( sugar ) out << "number of lists: ";
        out << vecolist.size();
        for( int subsim = 0; subsim < vecolist.size(); subsim++ ) {
            if( sugar ) out << subsim << ": ";
            out << vecolist.at(subsim).size() << " ";
            if( sugar ) out << "- ";
            for( const int& supsim : vecolist.at(subsim) )
                out << supsim << " ";
            // for( int supsim = 0; supsim < vecolist.at(subsim).size(); supsim++ ) {
                // out << vecolist.at(subsim).at(supsim) << " ";
            // }
            out << std::endl;
        }
    }
    
    
    /* Finished */
}

SimplicialMesh readSimplicialMeshPath( std::istream& in )
{
    assert( false );
    int innerdimension;
    int outerdimension;
    /* Preamble */
    in >> innerdimension >> outerdimension;
    /* Coordinates */
    Coordinates coords = readCoordinates( in );
    /* Combinatorial data */
    
    /* Combinatorial data - count of subsimplices */
    std::vector<int> simplexcounter(innerdimension+1);
    for( int d = 0; d <= innerdimension; d++ )
        in >> simplexcounter.at(d);
    
    /* Combinatorial data - subsimplices */
    std::map< std::pair<int,int>, std::vector<IndexMap> > sub;
    int subN;
    in >> subN;
    for( int n = 0; n < subN; n++ ) {
        
        std::pair<int,int> fromto;
        in >> fromto.first >> fromto.second;
        int numentries;
        in >> numentries;
        
        std::vector<IndexMap> vecoim( numentries, 
                                      IndexMap( IndexRange(0, countsubsimplices(fromto.first,fromto.second) ), 
                                                IndexRange(0, simplexcounter.at(fromto.second) )
                                              )
                                    );
        
        for( int supsim = 0; supsim < vecoim.size(); supsim++ ) {
            
            int CSS = countsubsimplices( fromto.first, fromto.second );
            
            for( int sub = 0; sub < CSS; sub++ )
                in >> vecoim.at(supsim).at( sub );
            
        }
        
        sub.insert( std::make_pair( fromto, vecoim ) );
        
    }
    
    /* Combinatorial data - supersimplices */
    std::map< std::pair<int,int>, std::vector<std::list<int>> > super;
    int superN;
    in >> superN;
    for( int n = 0; n < superN; n++ ) {
        std::pair<int,int> fromto;
        in >> fromto.first >> fromto.second;
        int numentries;
        in >> numentries;
        
        std::vector<std::list<int>> vecolist(numentries);
        
        for( int subsim = 0; subsim < vecolist.size(); subsim++ ) {
        
            auto& suplist = vecolist.at(subsim);
            int M;
            in >> M;
            for( int m = 0; m < M; m++ ) {
                int v;
                in >> v;
                suplist.push_back( v );
            }
            
        }
        
        super.insert( std::make_pair( fromto, vecolist ) );
        
    }
    
    /* Finish */
    return SimplicialMesh( innerdimension, outerdimension, coords, sub, super );
    
}



        
        
        
        
        