#ifndef INCLUDEGUARD_UTILITY_PROFILERUTILS_HPP
#define INCLUDEGUARD_UTILITY_PROFILERUTILS_HPP

#include <string>
#include <vector>

#include "../basic.hpp"

class SectionProfiler
{
    
    private:
    
        std::vector<  timestamp> times;
        std::vector<std::string> texts;
        
    public:
        
        explicit SectionProfiler( std::string text = "---" )
        {
            ping( text );
        }
        
        virtual ~SectionProfiler() {

            std::string text;

            int max_text_len = texts.at(0).size();
            for( int i = 1; i < times.size(); i++ )
                max_text_len = maximum( max_text_len, SIZECAST( texts.at(i).size() ) );

            for( int i = 1; i < times.size(); i++ ) {
                LOGPRINTF( "%12s %*s:\n",
                           timestamp2measurement( times.at(i) - times.at(i-1) ).c_str(),
                           max_text_len, texts[i].c_str()
                         );
            }
            LOGPRINTF("\n");
            
        }

        void ping( std::string text = "---" ) 
        {
            times.push_back( gettimestamp() );
            texts.push_back( text           );
        }
        
};



class StopWatch {

    private:

        timestamp   start_time;
        std::string text;

    public:

        StopWatch( std::string text = "---" ) {
            start_time = gettimestamp();
        }

        ~StopWatch() {
            timestamp end_time = gettimestamp();
            LOGPRINTF( "%12s %*s:\n",
                        timestamp2measurement( end_time - start_time ).c_str(),
                        static_cast<int>(text.length()), text.c_str()
                     );
        }
};


#endif
