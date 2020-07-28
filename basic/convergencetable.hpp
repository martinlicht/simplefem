#ifndef INCLUDEGUARD_BASIC_CONVERGENCETABLE
#define INCLUDEGUARD_BASIC_CONVERGENCETABLE

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    
    public:
        
        ConvergenceTable& operator<<( Float entry )
        {
            
            if( entries.size() == 0 )
                entries.push_back( std::vector<Float>(0) );
            
            entries.back().push_back( entry );
            
            return *this;
        }
        
        ConvergenceTable& operator<<( char code )
        {
            
            if( code == nl )
                entries.push_back( std::vector<Float>(0) );
            
            return *this; 
        }
        
        void print( std::ostream& out ) const 
        {
            
            auto temp = out.precision();
            out.precision(3);
            
            for( int i = 0; i < entries.size(); i++ )
            {
            
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    out << entries[i][j] << tab;
                    if( i == 0 )
                        out << "--";
                    else 
                        out << log2( entries[i-1][j] / entries[i][j] );
                    out << tab;
                }        
                
                out << nl;
                
            }
            
            out.precision(temp);
            
        }

    private:
        
        std::vector<std::vector<Float>> entries;
        
    
};




#endif
