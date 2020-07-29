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


inline std::string experimentfile( const std::string& basename )
{
    std::string ret;
    for( int i = 0; fileexists( ret = basename + std::string(".") + std::to_string(i) + std::string(".vtk") ); i++ );
    return ret;
}

inline std::string getbasename( const std::string& path )
{
      std::size_t begin = path.find_last_of("/");
      std::size_t   end = path.find_last_of(".");
      assert( begin != std::string::npos );
      assert( end   != std::string::npos );
      assert( end > begin+1 );
      if( begin == std::string::npos ) begin = 0;
      if( end   == std::string::npos ) end   = path.length();
      
      return path.substr(begin+1,end-begin-1);
}



#endif
