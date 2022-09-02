
#include <cstdio>
#include <ostream>
// #include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "../basic.hpp"

#include "convergencetable.hpp"



ConvergenceTable::ConvergenceTable( std::string table_name )
: table_name(table_name), 
    display_convergence_rates( true ),
    print_rowwise_instead_of_columnwise(false)
{
    make_new_row = true;
}
        
ConvergenceTable& ConvergenceTable::operator<<( Float entry )
{
    
    if( make_new_row ) {
        entries.push_back( std::vector<Float>(0) );
        make_new_row = false;
    }
    
    entries.back().push_back( entry );
    
    return *this;
}
        
ConvergenceTable& ConvergenceTable::operator<<( const std::string& seriesheader )
{
    
    seriesheaders.push_back( seriesheader );
    
    return *this;
}
        
        
ConvergenceTable& ConvergenceTable::operator<<( char code )
{
    
    if( code == nl )
        make_new_row = true; 
    else
        assert(false);
    
    return *this; 
}


std::string ConvergenceTable::text() const
{
    return text( display_convergence_rates );
}
        
void ConvergenceTable::lg() const
{
    LOG << text();
}

void ConvergenceTable::lg( bool display_convergence_rates ) const
{
    LOG << text( display_convergence_rates );
}

void ConvergenceTable::print( std::ostream& os )
{
    os << text();
}
        

std::string ConvergenceTable::text( bool display_convergence_rates ) const
{
    std::ostringstream ss;

    ss << std::string( 80, '-' ) << nl;
    
    if( print_rowwise_instead_of_columnwise )
        ss << text_transpose( display_convergence_rates );
    else
        ss << text_standard( display_convergence_rates );
    
    ss << nl;
    
    ss << std::string( 80, '-' ) << nl;

    return ss.str();
}

std::string ConvergenceTable::text_standard( bool display_convergence_rates ) const
{
    
    std::ostringstream ss;

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
    printf_into_stream( ss,  "\n%s\n", table_name.c_str() );

    if( entries.empty() ) {
        printf_into_stream( ss,  "----------- Table is empty!\n" );
        return ss.str();
    }
    
    const int num_series = entries[0].size();

    // if necessary, print column headers 
    if( not seriesheaders.empty() )
    {
        
        Assert( seriesheaders.size() == num_series, seriesheaders.size(), num_series );

        printf_into_stream( ss,  "%s %s", std::string( nc_indent_width, ' ' ).c_str(), column_separator ); // printf_into_stream( ss, "   \t");
            
        for( int j = 0; j < seriesheaders.size(); j++ )
        {
            
            std::string seriesheader = seriesheaders[j];
            
            if( seriesheader.size() > nc_cell_width ) {
                seriesheader.resize( nc_cell_width );
                seriesheader[ nc_cell_width-1 ] = '~';
            }
            
            printf_into_stream( ss,  "%*s%s", nc_cell_width, seriesheader.c_str(), column_separator );
            
            if( display_convergence_rates ) {
                printf_into_stream( ss,  "%s%s", std::string( nc_rate_width, ' ' ).c_str(), column_separator ); // printf_into_stream( ss, "          \t");
            }

        }
        
        printf_into_stream( ss, "\n");
        
    }
    
    // print the entries of the table, row by row 
    for( int i = 0; i < entries.size(); i++ )
    {
        
        printf_into_stream( ss,  "%*d:%s", nc_indent_width, i, column_separator );

        assert( entries[i].size() == num_series );
        
        // in the current row, print the entries 
        for( int j = 0; j < entries[i].size(); j++ )
        {
            
            printf_into_stream( ss, "% *.*e%s", nc_cell_width, nc_cell_precision, (double) entries[i][j], column_separator ); 
            
            if( display_convergence_rates ){
                
                if( i == 0 ) {
                    
                    printf_into_stream( ss,  "%s", std::string( nc_rate_width, '-' ).c_str() ); //printf_into_stream( ss, "----------");
                
                } else {
                
                    if( entries[i][j] > 0. and entries[i-1][j] > 0. ) {

                        double computed_rate = (double)std::log2( entries[i-1][j] / entries[i][j] );
                        
                        if( rates_are_float ) { 
                            printf_into_stream( ss, "%*.*e", nc_rate_width, nc_rate_precision, computed_rate  );
                        } else {
                            printf_into_stream( ss, "%*.*f", nc_rate_width, nc_rate_precision, computed_rate  );
                        }

                    } else {
                        
                        printf_into_stream( ss,  "%s", std::string( nc_rate_width, '$' ).c_str() ); //printf_into_stream( ss,  "%s", "$$$$$$$$$$" );

                    }
                
                }
                
                printf_into_stream( ss,  "%s", column_separator );
            }
            
        }        
        
        printf_into_stream( ss, "\n");
        
    }

    return ss.str();
                
}


std::string ConvergenceTable::text_transpose( bool display_convergence_rates ) const
{
    
    std::ostringstream ss;

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
    printf_into_stream( ss,  "\n%s\n", table_name.c_str() );

    const int num_entries_per_series = entries.size(); 
    
    if( entries.empty() ) {
        printf_into_stream( ss,  "----------- Table is empty!\n" );
        return ss.str();
    }
    
    const int num_series = entries[0].size();

    for( int j = 0; j < num_series; j++ )
    {
        // print the name of series first 
        if( not seriesheaders.empty() )
        {
        
            assert( seriesheaders.size() == num_series );

            std::string seriesheader = seriesheaders[j];

            if( seriesheader.size() > nc_header_width ) {
                seriesheader.resize( nc_header_width );
                seriesheader[ nc_header_width-1 ] = '~';
            }

            printf_into_stream( ss,  "%*s%s", nc_header_width, seriesheader.c_str(), cell_separator );

        }

        for( int i = 0; i < num_entries_per_series; i++ ) {

            assert( entries[i].size() == num_series );
            
            printf_into_stream( ss, "% *.*e%s", nc_cell_width, nc_cell_precision, (double) entries[i][j], cell_separator );
            
        }

        printf_into_stream( ss, "\n");
        
        if( display_convergence_rates )
        {
            
            if( not seriesheaders.empty() )
                printf_into_stream( ss,  "%s%s", std::string( nc_header_width, ' ' ).c_str(), cell_separator );

            for( int i = 0; i < num_entries_per_series; i++ )
            {

                if( i == 0 ) {
                        
                    printf_into_stream( ss,  "%s%s", std::string( nc_rate_width, '-' ).c_str(), cell_separator ); //printf_into_stream( ss, "----------");
                    
                } else {
                
                    if( entries[i][j] > 0. and entries[i-1][j] > 0. ) {

                        double computed_rate = std::log2( entries[i-1][j] / entries[i][j] );
                        
                        if( rates_are_float ) { 
                            printf_into_stream( ss, "% *.*e%s", nc_rate_width, nc_rate_precision, computed_rate, cell_separator );
                        } else {
                            printf_into_stream( ss, "% *.*f%s", nc_rate_width, nc_rate_precision, computed_rate, cell_separator );
                        }

                    } else {
                        
                        printf_into_stream( ss,  "%s%s", std::string( nc_rate_width, '$' ).c_str(), cell_separator ); //printf_into_stream( ss,  "%s", "$$$$$$$$$$" );

                    }
                
                }

            }

            printf_into_stream( ss, "\n");

        }

    }

    return ss.str();
            
}





//         void print_stream( std::ostream& ss, bool display_convergence_rates ) const 
//         {
//             std::ostringstream str;
//             
//             for( int i = 0; i < entries.size(); i++ )
//             {
//                 
//                 ss << std::setw(3) << i << ":\t";
// 
//                 assert( entries[i].size() == entries.front().size() );
//                 
//                 for( int j = 0; j < entries[i].size(); j++ )
//                 {
//                     
//                     ss << std::setprecision(6) << std::scientific << std::showpos << (double) entries[i][j];
//                     
//                     if( display_convergence_rates ){
//                         
//                         if( i == 0 ) {
//                             
//                             ss << "----------";
//                         
//                         } else {
//                         
//                             if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
//                                 ss << std::setw(10) << std::setprecision(3) << std::scientific << std::showpos << (double) std::log2( entries[i-1][j] / entries[i][j] );
//                             else
//                                 ss << "$$$$$$$$$$";
//                         
//                         }
//                         
//                         ss << '\t';
//                     }
//                     
//                 }        
//                 
//                 ss << nl;
//                 
//             }
//             
//             ss << str.str();
//             
//         }

        // Float get_convergence_rate( int row, int column )
        // {
        //     assert( 0 <= row );
        //     assert( row < entries.size() );
        //     assert( 1 <= row );
        //     assert( 0 <= column );
        //     assert( column < entries[row].size() );
        //     assert( column < entries[row-1].size() );
        //     Float ret = std::log2( entries[row-1][column] / entries[row][column] );
        //     return ret;
        // }




