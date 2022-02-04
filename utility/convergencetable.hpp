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

        explicit ConvergenceTable( std::string table_name = "---------- Default Table Name ----------", bool display_convergence_rates = true )
        : make_new_row(true), 
          table_name(table_name), 
          display_convergence_rates( display_convergence_rates ),
          print_transpose_instead_of_standard(false)
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



        void print( std::ostream& os, bool display_convergence_rates ) const
        {
            
            if( print_transpose_instead_of_standard )
                print_transpose( os, display_convergence_rates );
            else
                print_standard( os, display_convergence_rates );
        }
        
        
        
        // TODO 
        // Introduced temporarily until format library is available
        // C++ streams are currently not supported, 
        // instead use printf from the C library
        void print_standard( std::ostream&, bool display_convergence_rates ) const
        {
            
            
            // introduce several constants that drive the output format 
            
            const char* column_separator = "   ";

            const bool rates_are_float = false;
            
            const int nc_indent_width = 3;
            
            const int nc_cell_precision = 3;

            const int nc_rate_precision = 2;

            const int nc_cell_width = 7 + nc_cell_precision + 0;
            // 12; sign + digit + . + e + sign + two digits = 6 chars 

            const int nc_rate_width = ( rates_are_float ? ( 7 + nc_rate_precision ) : 3 + nc_rate_precision ) + 0;
            // if float: 10; sign + digit + . + e + sign + two digits = 7 chars 
            // if fixed:  6; sign + digit + . = 3 chars 

        
            // First line is the name of the table 
            std::printf( "\n%s\n", table_name.c_str() );

            if( entries.empty() ) {
                std::printf( "----------- Table is empty!\n" );
                return;
            }
            
            // if necessary, print column headers 
            if( not columnheaders.empty() )
            {
                
                std::printf( "%s %s", std::string( nc_indent_width, ' ' ).c_str(), column_separator ); // std::printf("   \t");
                    
                for( int j = 0; j < columnheaders.size(); j++ )
                {
                    
                    std::string columnheader = columnheaders[j];
                    
                    if( columnheader.size() > nc_cell_width ) {
                        columnheader.resize( nc_cell_width );
                        columnheader[ nc_cell_width-1 ] = '~';
                    }
                    
                    std::printf( "%*s%s", nc_cell_width, columnheader.c_str(), column_separator );
                    
                    if( display_convergence_rates ) {
                        std::printf( "%s%s", std::string( nc_rate_width, ' ' ).c_str(), column_separator ); // std::printf("          \t");
                    }

                }
                
                std::printf("\n");
                
            }
            
            // print the entries of the table, row by row 
            for( int i = 0; i < entries.size(); i++ )
            {
                
                std::printf( "%*d:%s", nc_indent_width, i, column_separator );

                assert( entries[i].size() == entries.front().size() );
                
                // in the current row, print the entries 
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    
                    std::printf("% *.*Le%s", nc_cell_width, nc_cell_precision, (long double) entries[i][j], column_separator ); 
                    
                    if( display_convergence_rates ){
                        
                        if( i == 0 ) {
                            
                            std::printf( "%s", std::string( nc_rate_width, '-' ).c_str() ); //std::printf("----------");
                        
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) {

                                long double computed_rate = std::log2( entries[i-1][j] / entries[i][j] );
                                
                                if( rates_are_float ) { 
                                    std::printf("%*.*Le", nc_rate_width, nc_rate_precision, computed_rate  );
                                } else {
                                    std::printf("%*.*Lf", nc_rate_width, nc_rate_precision, computed_rate  );
                                }

                            } else {
                                
                                std::printf( "%s", std::string( nc_rate_width, '$' ).c_str() ); //std::printf( "%s", "$$$$$$$$$$" );

                            }
                        
                        }
                        
                        std::printf( "%s", column_separator );
                    }
                    
                }        
                
                std::printf("\n");
                
            }
                        
        }


        // TODO 
        // Introduced temporarily until format library is available
        // C++ streams are currently not supported, 
        // instead use printf from the C library
        void print_transpose( std::ostream&, bool display_convergence_rates ) const
        {
            
            // introduce several constants that drive the output format 
            
            const char* cell_separator = "   ";

            const bool rates_are_float = false;
            
            const int nc_header_width = 10;
            
            const int nc_cell_precision = 3;

            const int nc_rate_precision = 2;

            const int nc_cell_width = 6 + nc_cell_precision + 2;
            // 12; sign + digit + . + e + sign + two digits = 6 chars 

            const int nc_rate_width = std::max( nc_cell_width, ( rates_are_float ? ( 7 + nc_rate_precision ) : 3 + nc_rate_precision ) + 0 );
            // see above ....
        
            // First line is the name of the table 
            std::printf( "\n%s\n", table_name.c_str() );

            const int num_entries_per_series = entries.size(); 
            
            if( entries.empty() ) {
                std::printf( "----------- Table is empty!\n" );
                return;
            }
            
            const int num_series = entries[0].size();

            for( int j = 0; j < num_series; j++ )
            {
                // print the name of series first 
                if( not columnheaders.empty() )
                {
                
                    assert( columnheaders.size() == num_series );

                    std::string columnheader = columnheaders[j];

                    if( columnheader.size() > nc_header_width ) {
                        columnheader.resize( nc_header_width );
                        columnheader[ nc_header_width-1 ] = '~';
                    }

                    std::printf( "%*s%s", nc_header_width, columnheader.c_str(), cell_separator );

                }

                for( int i = 0; i < num_entries_per_series; i++ ) {

                    assert( entries[i].size() == num_series );
                    
                    std::printf("% *.*Le%s", nc_cell_width, nc_cell_precision, (long double) entries[i][j], cell_separator );
                    
                }

                std::printf("\n");
                
                if( display_convergence_rates )
                {
                    
                    if( not columnheaders.empty() )
                        std::printf( "%s%s", std::string( nc_header_width, ' ' ).c_str(), cell_separator );

                    for( int i = 0; i < num_entries_per_series; i++ )
                    {

                        if( i == 0 ) {
                                
                            std::printf( "%s%s", std::string( nc_rate_width, '-' ).c_str(), cell_separator ); //std::printf("----------");
                            
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) {

                                long double computed_rate = std::log2( entries[i-1][j] / entries[i][j] );
                                
                                if( rates_are_float ) { 
                                    std::printf("% *.*Le%s", nc_rate_width, nc_rate_precision, computed_rate, cell_separator );
                                } else {
                                    std::printf("% *.*Lf%s", nc_rate_width, nc_rate_precision, computed_rate, cell_separator );
                                }

                            } else {
                                
                                std::printf( "%s%s", std::string( nc_rate_width, '$' ).c_str(), cell_separator ); //std::printf( "%s", "$$$$$$$$$$" );

                            }
                        
                        }

                    }

                    std::printf("\n");

                }

            }
                    
        }

    private:
        
        std::vector<std::vector<Float>> entries;
        std::vector<std::string> columnheaders;
        bool make_new_row;
        
    public:
        
        std::string table_name;
        bool display_convergence_rates;
        bool print_transpose_instead_of_standard;
    
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
