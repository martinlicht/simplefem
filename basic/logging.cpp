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
    
    //FILE* f = stderr; //use_cerr ? stderr : stdout;
    FILE* f = ( use_cerr ? stderr : stdout );
    
    // const auto str = this->str();
    const auto& str = internal;
    
    bool use_prefix_next  = log_has_a_fresh_line;
    
    if( str.empty() ) {
        // internalstream << prefix;
        // std::cout << "\nEMPTY\n";
        return;
    }

    std::string prefix = digitalcodenow(); 
    
    for( int c = 0; c < str.size(); c++ )
    {
        // "\033[46m %3d\033[m"
        if( use_prefix_next ) { 
            use_prefix_next = false;
            #ifdef USE_COLORED_OUTPUT
            fprintf( f, "\033[96m[%s %s %4u]\033[m\t", prefix.c_str(), filename.c_str(), linenumber );
            #else 
            fprintf( f, "[%s %s %4u]\t", prefix.c_str(), filename.c_str(), linenumber );
            #endif 
        }
        
        auto character = str.at(c);
        
        fputc( character, f ); // internalstream << character;

        if( character == '\n' ) 
        { 
            // internalstream.flush();
            use_prefix_next = true;
        }
        
    }
    
    log_has_a_fresh_line = false;
    
    if( not str.empty() && str.back() == '\n' ) {
        log_has_a_fresh_line = true;
    }
    
    if( not str.empty() && str.back() != '\n' && pad_newline_if_there_is_none ) {
        log_has_a_fresh_line = true;
        fputc( nl, f ); // internalstream << nl;                
    }
    
    

    if( print_file_and_line and log_has_a_fresh_line ) {
        fputs( prefix.c_str(), f ); // internalstream << prefix;
        // internalstream << "\e[91m" << filename << ':' << linenumber << "\e[39m" << '\n';
        fprintf( f, "%s:%d\n", filename.c_str(), linenumber ); // internalstream << "" << filename << ':' << linenumber << '\n';
    }

    //#ifndef NDEBUG
    // internalstream.flush();
    fflush(f);
    //#endif
    
}






#if defined(_OPENMP)
#include <omp.h>

OpenMP_Reporter::OpenMP_Reporter()
{
    
    #if defined(_OPENMP)
    LOG << "OpenMP Reporter: started\n";
    LOG << "\tOpenMP Value: " << _OPENMP << nl;
    LOG << "\tMaximum number of threads: " << omp_get_max_threads() << nl;
    LOG << "\tThread limit: " << omp_get_thread_limit() << nl;
    // LOG << "\tMaximum number of processors: " << p << " - > " << omp_get_place_num_procs() << nl;
    LOG << "\tMaximum number of places: " << omp_get_num_places() << nl;
    for( int p = 0; p < omp_get_num_places(); p++ ) {
        LOG << "\t\tMaximum number of processors per place: " << p << " - > " << omp_get_place_num_procs(p) << nl;
    }
    #endif
}

OpenMP_Reporter::~OpenMP_Reporter()
{
    #if defined(_OPENMP)
    LOG << "OpenMP Reporter: finished\n";
    #endif
}

OpenMP_Reporter  omp_reporter;

#endif // #if defined(_OPENMP)



#endif
