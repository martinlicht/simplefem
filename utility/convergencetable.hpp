#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

#include <cstdio>
#include <iostream>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    
    public:
        
        explicit ConvergenceTable()
        : make_new_row(true)
        {
            
        }
        
        ConvergenceTable& operator<<( Float entry )
        {
            
            if( make_new_row ) {
                entries.push_back( std::vector<Float>(0) );
                make_new_row = false;
            }
            
            entries.back().push_back( entry );
            
            return *this;
        }
        
        ConvergenceTable& operator<<( char code )
        {
            
            if( code == nl )
                make_new_row = true; //entries.push_back( std::vector<Float>(0) );
            
            return *this; 
        }
        
//         void print( std::ostream& out ) const 
//         {
//             
//             auto temp = out.precision(); // TODO include fmt as soon as available 
//             out.precision(myprecision);
//             
//             for( int i = 0; i < entries.size(); i++ )
//             {
//                 
//                 out << i << ":" << tab;
//                 
//                 for( int j = 0; j < entries[i].size(); j++ )
//                 {
//                     out << entries[i][j] << tab;
//                     if( i == 0 )
//                         out << "--";
//                     else 
//                         out << std::log2( entries[i-1][j] / entries[i][j] );
//                     out << tab;
//                 }        
//                 
//                 out << nl;
//                 
//             }
//             
//             out.precision(temp);
//             
//         }
        
        void print( std::ostream&, bool show_rates = true ) const // TODO introduced temporarily until format library is available
        {
            
            for( int i = 0; i < entries.size(); i++ )
            {
                
                std::printf("%3d:\t",i);

                assert( entries[i].size() == entries.front().size() );
                
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    
                    std::printf("%.6Le\t", (long double) entries[i][j] ); 
                    
                    if( show_rates ){
                        
                        if( i == 0 ) {
                            
                            std::printf("        --");
                        
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
                                std::printf("%10.3Le", (long double) std::log2( entries[i-1][j] / entries[i][j] ) );
                            else
                                std::printf("        $$");
                        
                        }
                        
                        std::printf("\t");
                    }
                    
                }        
                
                std::printf("\n");
                
            }
                        
        }

    private:
        
        std::vector<std::vector<Float>> entries;
        
        bool make_new_row;
    
};




#endif
