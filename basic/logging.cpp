#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING
bool log_has_a_fresh_line = true;

// #include <iostream>
// #include <iomanip>
#include <cstdio>

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

        

Logger::~Logger()
{
    
    FILE* f = use_cerr ? stderr : stdout;
    
    const auto& str = *this; //this->str();
    
    bool use_prefix_next  = log_has_a_fresh_line;
    
    if( str.empty() ) {
        // internalstream << prefix;
        // std::cout << "\nEMPTY\n";
        return;
    }

    std::string prefix = protocolprefixnow();
    
    for( int c = 0; c < str.size(); c++ )
    {

        if( use_prefix_next ) { 
            use_prefix_next = false;
            fputs( prefix.c_str(), f ); // internalstream << prefix;
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
