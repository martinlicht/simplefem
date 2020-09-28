#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

#include <cstdio>
#include <iostream>
#include <iomanip>
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
        
        ConvergenceTable& operator<<( char code )
        {
            
            if( code == nl )
                make_new_row = true; //entries.push_back( std::vector<Float>(0) );
            
            return *this; 
        }
        
        void print( std::ostream& os )
        {
            print( os, display_convergence_rates );
        }
        
        void print( std::ostream&, bool display_convergence_rates ) const // TODO introduced temporarily until format library is available
        {
            
            for( int i = 0; i < entries.size(); i++ )
            {
                
                std::printf("%3d:\t",i);

                assert( entries[i].size() == entries.front().size() );
                
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    
                    std::printf("%.6Le\t", (long double) entries[i][j] ); 
                    
                    if( display_convergence_rates ){
                        
                        if( i == 0 ) {
                            
                            std::printf("----------");
                        
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
                                std::printf("%10.3Le", (long double) std::log2( entries[i-1][j] / entries[i][j] ) );
                            else
                                std::printf("$$$$$$$$$$");
                        
                        }
                        
                        std::printf("\t");
                    }
                    
                }        
                
                std::printf("\n");
                
            }
                        
        }


        void print_iostream( std::ostream& os, bool display_convergence_rates ) const 
        {
             // TODO use {fmt} library as soon as available 
            std::ostringstream str;
            
            for( int i = 0; i < entries.size(); i++ )
            {
                
                os << std::setw(3) << i << ":\t";

                assert( entries[i].size() == entries.front().size() );
                
                for( int j = 0; j < entries[i].size(); j++ )
                {
                    
                    os << std::setprecision(6) << std::scientific << std::showpos << (long double) entries[i][j];
                    
                    if( display_convergence_rates ){
                        
                        if( i == 0 ) {
                            
                            os << "----------";
                        
                        } else {
                        
                            if( entries[i][j] > 0. and entries[i-1][j] > 0. ) 
                                os << std::setw(10) << std::setprecision(3) << std::scientific << std::showpos << (long double) std::log2( entries[i-1][j] / entries[i][j] );
                            else
                                os << "$$$$$$$$$$";
                        
                        }
                        
                        os << '\t';
                    }
                    
                }        
                
                os << '\n';
                
            }
            
            os << str.str();
            
        }


    private:
        
        std::vector<std::vector<Float>> entries;
        
        bool make_new_row;
        bool display_convergence_rates;
    
};




#endif
