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

inline bool fileexists( std::string filename )
{
    std::ifstream file( filename.c_str() );
    return file.good();
}

inline std::string adaptfilename( std::string filename )
{
    if( !fileexists( filename ) ) 
        return filename;
    
    auto dotpos = filename.rfind( "." );
    std::string pre_dot  = filename.substr(0,dotpos);
    std::string post_dot = filename.substr(dotpos+1);
    
    for( int i = 0; fileexists( filename = pre_dot + std::string(".") + std::to_string(i) + std::string(".") + post_dot ); i++ );
    
    return filename;
}

  


#endif
