#ifndef INCLUDEGUARD_UTILITY_TIMEBEACON_HPP
#define INCLUDEGUARD_UTILITY_TIMEBEACON_HPP

#include <string>
#include <vector>

#include "../basic.hpp"

class TimeBeacon
{
    
    private:
    
        std::vector<  timestamp> times;
        std::vector<std::string> texts;
        
    public:
        
        explicit TimeBeacon( std::string text = "---" )
        {
            ping( text );
        }
        
        void ping( std::string text = "---" ) 
        {
            times.push_back( gettimestamp() );
            texts.push_back( text           );
        }
        
        void show() const {
            
            std::string text;
//             for( int i = 1; i < times.size(); i++ )
//                 text += texts.at(i) + " -> " + timestamp2measurement( times.at(i) - times.at(i-1) ) +  "\t";
//             LOG << text << nl;
            int max_text_len = texts.at(0).size();
            for( int i = 1; i < times.size(); i++ ) max_text_len = maximum( max_text_len, (int)texts.at(i).size() );
            for( int i = 1; i < times.size(); i++ )
                LOGPRINTF( "%*s: %12s\t",
                           max_text_len, texts[i].c_str(),
                           timestamp2measurement( times.at(i) - times.at(i-1) ).c_str()
                         );
//                 LOG << " ----> " << timestamp2measurement( times.at(i) - times.at(i-1) ) << " @ " << texts.at(i) << nl;
                LOGPRINTF("\n");
            
        }
        
        virtual ~TimeBeacon() {
            show();
        }

};

#endif
