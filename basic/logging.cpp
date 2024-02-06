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
    
    
    if( str.empty() ) return;

    
    std::string prefix = digitalcodenow(); 
    
    for( int c = 0; c < str.size(); c++ )
    {
    
        if( use_prefix_next ) { 
    
            use_prefix_next = false;
    
            std::string formatstring("[%s %s %4u]\t");
            #ifdef USE_COLORED_OUTPUT
            formatstring = "\033[96m" + formatstring + "\033[m";
            #endif 
            fprintf( f, formatstring.c_str(), prefix.c_str(), filename.c_str(), linenumber );
    
        }
        
        auto character = str.at(c);
        fputc( character, f );

        if( character == '\n' ) 
            use_prefix_next = true;
        
    }
    
    log_has_a_fresh_line = false;
    
    if( not str.empty() && str.back() == '\n' ) {
        log_has_a_fresh_line = true;
    }
    
    if( not str.empty() && str.back() != '\n' && pad_newline_if_there_is_none ) {
        log_has_a_fresh_line = true;
        fputc( nl, f ); 
    }
    
    

    if( print_file_and_line and log_has_a_fresh_line ) {
        fputs( prefix.c_str(), f ); 
        fprintf( f, "%s:%d\n", filename.c_str(), linenumber ); 
    }

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
