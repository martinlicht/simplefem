#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING
bool log_has_a_fresh_line = true;

// #include <iostream>
// #include <iomanip>
#include <cstdio>



enum class TextColors
{
    black          = 30,
    black_light    = 90,
    
    red            = 31,
    red_light      = 91,
    
    green          = 32,
    green_light    = 92, 
    
    yellow         = 33, 
    yellow_bright  = 93, 
    
    blue           = 34,
    blue_bright    = 94,
    
    magenta        = 35,
    magenta_bright = 95,
    
    cyan           = 36,
    cyan_light     = 96, 
    
    white          = 37,
    white_light    = 97 
};


Logger::Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline,
            const char* filename,
            const int linenumber
        )
        : 
        use_cerr( use_cerr ), // internalstream( use_cerr ? std::cerr : std::cout ), // internalstream( os ),
        pad_newline_if_there_is_none( do_newline ),
        filename( filename ),
        linenumber( linenumber )
        {
            //*this << std::setprecision(10);
        }


auto digitalcodenow() -> std::string;

Logger::~Logger()
{
    
    FILE* f = ( use_cerr ? stderr : stdout );
    
    const auto& str = internal;
    
    bool use_prefix_next  = log_has_a_fresh_line;
    
    
    // if the log string is empty, just do not print anything 

    if( str.empty() ) return;

    
    // complete the prefix string 

    const std::string time_code = digitalcodenow(); 

    const std::string formatstring("[%s %s %4u]\t");
    #ifdef USE_COLORED_OUTPUT
    const std::string colorcode_begin = "\033[96m";
    const std::string colorcode_close = "\033[m";
    #else
    const std::string colorcode_begin = "";
    const std::string colorcode_close = "";
    #endif
    const std::string prefix = printf_into_string(
        "%s[%s %s %4d]%s\t", 
        colorcode_begin.c_str(), time_code.c_str(), filename.c_str(), linenumber, colorcode_close.c_str() );
    
    
    // assemble final output string 

    int number_of_newlines = 0;
    for( char character : str ) 
        if( character == '\n' )
            number_of_newlines++;

    std::string output_string;

    output_string.reserve( 1 + str.size() + number_of_newlines * prefix.size() );

    for( char character : str ) {
    
        if( use_prefix_next ) { 
    
            use_prefix_next = false;
    
            output_string += prefix;
            // fputs( prefix.c_str(), f );
            
        }
        
        output_string += character;
        // fputc( character, f );

        if( character == '\n' ) 
            use_prefix_next = true;
        
    }
    
    log_has_a_fresh_line = false;
    
    if( str.back() == '\n' or pad_newline_if_there_is_none ) {
        
        log_has_a_fresh_line = true;
        
        if( str.back() != '\n' and pad_newline_if_there_is_none ) 
            output_string += nl;
            // fputc( nl, f ); 
    
    }

    fputs( output_string.c_str(), f );
    fflush(f);
    
}






#if defined(_OPENMP)
#include <omp.h>

OpenMP_Reporter::OpenMP_Reporter()
{
    
    #if defined(_OPENMP)
    LOG << "###OMP###\tOpenMP Reporter: started\n";
    LOG << "###OMP###\tOpenMP Value: " << _OPENMP << nl;
    LOG << "###OMP###\tMaximum number of threads: " << omp_get_max_threads() << nl;
    LOG << "###OMP###\tThread limit: " << omp_get_thread_limit() << nl;
    // LOG << "###OMP###\tMaximum number of processors: " << p << " -> " << omp_get_place_num_procs() << nl;
    LOG << "###OMP###\tMaximum number of places: " << omp_get_num_places() << nl;
    for( int p = 0; p < omp_get_num_places(); p++ ) {
        LOG << "###OMP###\t\tMaximum number of processors per place: " << p << " - > " << omp_get_place_num_procs(p) << nl;
    }
    #endif
}

OpenMP_Reporter::~OpenMP_Reporter()
{
    #if defined(_OPENMP)
    LOG << "###OMP###\tOpenMP Reporter: finished\n";
    #endif
}

OpenMP_Reporter  omp_reporter;

#endif // #if defined(_OPENMP)



#endif
