#ifndef INCLUDEGUARD_UTILITY_FILES_HPP
#define INCLUDEGUARD_UTILITY_FILES_HPP


#include <fstream>
#include <string>




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


inline std::string experimentfile( const std::string& basename, const std::string extension = "vtk" )
{
    std::string ret;
    for( int i = 0; fileexists( ret = basename + std::string(".") + std::to_string(i) + std::string(".") + std::string(extension) ); i++ );
    return ret;
}

inline std::string getbasename( const std::string& path )
{
      std::string::size_type last_slash = path.find_last_of("/");
      std::string::size_type last_dot   = path.find_last_of(".");
      
      assert( last_dot != std::string::npos );
      
      std::string::size_type begin = ( last_slash != std::string::npos ) ? ( last_slash+1 ) : 0;
      std::string::size_type end   = last_dot;
      
      return path.substr(begin,end-begin);
}



#endif
