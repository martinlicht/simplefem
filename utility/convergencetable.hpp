#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    
    public:

        explicit ConvergenceTable( bool display_convergence_rates = true )
        : make_new_row(true), 
          display_convergence_rates( display_convergence_rates )
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
        
        ConvergenceTable& operator<<( const std::string& columnheader )
        {
            
            columnheaders.push_back( columnheader );
            
            return *this;
        }
        
        
        ConvergenceTable& operator<<( char code )
        {
            
            if( code == nl )
                make_new_row = true; 
            else
                assert(false);
            
            return *this; 
        }
        
        std::string text() const
        {
            return text( display_convergence_rates );
        }
        
        std::string text( bool display_convergence_rates ) const
        {
            std::ostringstream ss;
            print( ss, display_convergence_rates );
            return ss.str();
        }
        
        void lg() const
        {
            lg( display_convergence_rates );
        }
        
        void lg( bool display_convergence_rates ) const
        {
            LOG << text( display_convergence_rates );
        }
        
        void print( std::ostream& os )
        {
            print( os, display_convergence_rates );
        }
        


        Float get_convergence_rate( int row, int column )
        {
            assert( 0 <= row );
            assert( row < entries.size() );
            assert( 1 <= row );
            assert( 0 <= column );
            assert( column < entries[row].size() );
            assert( column < entries[row-1].size() );
            Float ret = std::log2( entries[row-1][column] / entries[row][column] );
            return ret;
        }

        
        
        const int tablecell_width = 12;

        const int tablecell_precision = 6;

        const int convergence_rate_width = 10;

        const int convergence_rate_precision = 3;

        const int initial_tab = 3;
        
        // TODO 
        // Introduced temporarily until format library is available
        // C++ streams are currently not supported, 
        // instead use printf from the C library
        void print( std::ostream&, bool display_convergence_rates ) const
        {
            
            // if necessary, print column headers 
            if( not columnheaders.empty() )
            {
                
                std::printf( "%s\t", std::string( 3, ' ' ).c_str() ); // std::printf("   \t");
                    
                for( int j = 0; j < columnheaders.size(); j++ )
                {
                    
                    std::string columnheader = columnheaders[j];
                    
                    if( columnheader.size() > 12 ) {
                        columnheader.resize(12);
                        columnheader[12] = '~';
                    }
                    
                    std::printf( "%*s\t", 12, columnheader.c_str() );
                    
                    if( display_convergence_rates ) 
                        std::printf( "%s\t", std::string(10,' ').c_str() ); // std::printf("          \t");
                    
                }
                
                std::printf("\n");
                
            }
            
            // print the entries of the table, row by row 
            for( int i = 0; i < entries.size(); i++ )
            {
                
                std::printf( "%*d:\t", 3, i );

                assert( entries[i].size() == entries.front().size() );
                
                // in the current row, print the entries 
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    
                    std::printf("%*.*Le\t", 12, 6, (long double) entries[i][j] ); 
                    
                    if( display_convergence_rates ){
                        
                        if( i == 0 ) {
                            
                            std::printf( "%s", std::string( 10, '-' ).c_str() ); //std::printf("----------");
                        
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
                                std::printf("%*.*Le", 10, 3, (long double) std::log2( entries[i-1][j] / entries[i][j] ) );
                            else
                                std::printf( "%s", std::string( 10, '$' ).c_str() ); //std::printf( "%s", "$$$$$$$$$$" );
                        
                        }
                        
                        std::printf("\t");
                    }
                    
                }        
                
                std::printf("\n");
                
            }
                        
        }

    private:
        
        std::vector<std::vector<Float>> entries;
        std::vector<std::string> columnheaders;
        
        bool make_new_row;
        bool display_convergence_rates;
    
};



//         void print_stream( std::ostream& os, bool display_convergence_rates ) const 
//         {
//              TODO use {fmt} library as soon as available 
//             std::ostringstream str;
//             
//             for( int i = 0; i < entries.size(); i++ )
//             {
//                 
//                 os << std::setw(3) << i << ":\t";
// 
//                 assert( entries[i].size() == entries.front().size() );
//                 
//                 for( int j = 0; j < entries[i].size(); j++ )
//                 {
//                     
//                     os << std::setprecision(6) << std::scientific << std::showpos << (long double) entries[i][j];
//                     
//                     if( display_convergence_rates ){
//                         
//                         if( i == 0 ) {
//                             
//                             os << "----------";
//                         
//                         } else {
//                         
//                             if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
//                                 os << std::setw(10) << std::setprecision(3) << std::scientific << std::showpos << (long double) std::log2( entries[i-1][j] / entries[i][j] );
//                             else
//                                 os << "$$$$$$$$$$";
//                         
//                         }
//                         
//                         os << '\t';
//                     }
//                     
//                 }        
//                 
//                 os << nl;
//                 
//             }
//             
//             os << str.str();
//             
//         }

#endif
