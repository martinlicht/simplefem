#ifndef INCLUDEGUARD_BASIC_CONVERGENCETABLE
#define INCLUDEGUARD_BASIC_CONVERGENCETABLE

#include <cstdio>
#include <iostream>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    
    public:
        
        explicit ConvergenceTable( unsigned int prec = 7 )
        : make_new_row(true), myprecision(prec)
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
        
        void print( std::ostream& ) const // TODO introduced temporarily until format library is available
        {
            
            for( int i = 0; i < entries.size(); i++ )
            {
                
                std::printf("%3d:\t",i);
                
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    std::printf("%.6Le\t", entries[i][j] ); 
                    if( i == 0 )
                        std::printf("        --");
                    else 
                        std::printf("%10.3Le", std::log2( entries[i-1][j] / entries[i][j] ) );
                    std::printf("\t");
                }        
                
                std::printf("\n");
                
            }
                        
        }

    public:
        
        std::vector<std::vector<Float>> entries;
        
        bool make_new_row;
        
        unsigned int myprecision;
    
};




#endif
