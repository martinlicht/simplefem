#ifndef INCLUDEGUARD_UTILITY_PROFILERUTILS_HPP
#define INCLUDEGUARD_UTILITY_PROFILERUTILS_HPP

#include <chrono>
#include <string>
#include <vector>

#include "../basic.hpp"

class SectionProfiler
{
    private:
    
        typedef std::chrono::steady_clock::time_point time_type;
    
        std::vector<  time_type> times;
        std::vector<std::string> texts;
        
    public:
        
        explicit SectionProfiler( std::string text = "---" )
        {
            ping( text );
        }
        
        virtual ~SectionProfiler() {

            ping("FINISH");

            std::string text;

            assert( texts.size() == times.size() and texts.size() >= 2 );

            for( int i = 0; i < times.size()-1; i++ ) {
                
                auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>( times.at(i+1) - times.at(i) ).count();

                LOGPRINTF( "%ju ns \t %s\n", static_cast<uintmax_t>( elapsed_time ), texts[i].c_str() );
            }
            LOGPRINTF("\n");
            
        }

        void ping( std::string text = "---" ) 
        {
            times.push_back( std::chrono::steady_clock::now() );
            texts.push_back( text                                      );
        }
        
};



class StopWatch {

    private:

        typedef std::chrono::steady_clock::time_point time_type;
    
        time_type   start_time;
        std::string text;

    public:

        StopWatch( std::string text = "---" ) 
        : start_time(std::chrono::steady_clock::now()), text(text) 
        {}

        ~StopWatch() {
            
            time_type end_time = std::chrono::steady_clock::now();
            
            auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time ).count();

            LOGPRINTF( "%ju ns \t %s\n", static_cast<uintmax_t>( elapsed_time ), text.c_str() );
            
        }
};


#endif
