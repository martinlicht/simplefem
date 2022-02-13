#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING
bool log_has_a_fresh_line = true;

#include <iostream>
#include <iomanip>

Logger::Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline,
            const char* filename,
            const int linenumber
        )
        : 
//         internalstream( os ),
        internalstream( use_cerr ? std::cerr : std::cout ),
        pad_newline_if_there_is_none( do_newline ),
        filename( filename ),
        linenumber( linenumber )
        {
            *this << std::setprecision(10);
        }

        

Logger::~Logger()
{
    
    const auto str = this->str();
    
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
            internalstream << prefix;
        }
        
        auto character = str.at(c);
        
        internalstream << character;

        if( character == '\n' ) 
        { 
            internalstream.flush();
            use_prefix_next = true;
        }
        
    }
    
    log_has_a_fresh_line = false;
    
    if( not str.empty() && str.back() == '\n' ) {
        log_has_a_fresh_line = true;
    }
    
    if( not str.empty() && str.back() != '\n' && pad_newline_if_there_is_none ) {
        log_has_a_fresh_line = true;
        internalstream << nl;                
    }
    
    

    if( print_file_and_line and log_has_a_fresh_line ) {
        internalstream << prefix;
        // internalstream << "\e[91m" << filename << ':' << linenumber << "\e[39m" << '\n';
        internalstream << "" << filename << ':' << linenumber << '\n';
    }

    #ifndef NDEBUG
    internalstream.flush();
    #endif
    
}

#endif
