

#include <algorithm>
#include <fstream>
#include <istream>
#include <ostream>
#include <map>
#include <string>
#include <utility>
#include <vector>



#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/io.floatvector.hpp"




void writeFloatVector( const char* filename, const FloatVector& vec )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeFloatVector( myfile, vec );
    myfile.close();
}

FloatVector readFloatVector( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    FloatVector vec = readFloatVector( myfile );
    myfile.close();
    return vec;
}



void writeFloatVector( std::ostream& out, const FloatVector& vec )
{
    out << vec.getdimension() << std::endl;
    for( int p = 0; p < vec.getdimension(); p++ )
      out << vec.at(p) << std::endl;
}

FloatVector readFloatVector( std::istream& in )
{
    int dim;
    in >> dim;
    FloatVector vec( dim );
    for( int p = 0; p < vec.getdimension(); p++ )
      in >> vec.at(p);
    return vec;
}


