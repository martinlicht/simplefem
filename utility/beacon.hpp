#ifndef INCLUDEGUARD_UTILITY_TIMEBEACON_HPP
#define INCLUDEGUARD_UTILITY_TIMEBEACON_HPP

#include <iostream>
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
            
            std::cout << ">>> " << timestamp2measurement( times.at(0) ) << " @ " << texts.at(0) << std::endl;
            for( int i = 1; i < times.size(); i++ )
                std::cout << ">>> " << timestamp2measurement( times.at(i) - times.at(i-1) ) << " @ " << texts.at(i) << std::endl;
            
        }
        
        virtual ~TimeBeacon() {
            show();
        }

        
    
};




#endif
