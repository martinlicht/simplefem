#ifndef INCLUDEGUARD_FILES_HPP
#define INCLUDEGUARD_FILES_HPP


#include <fstream>




inline std::fstream openinputfile( std::string filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    return myfile;
}

inline std::fstream openoutputfile( std::string filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    return myfile;
}


  


#endif
