#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

#include <ostream>
#include <string>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    
    private:
        
        std::vector<std::vector<Float>> entries;
        std::vector<std::string> seriesheaders;
        bool make_new_row;
        
    public:
        
        std::string table_name;
        bool display_convergence_rates;
        bool print_rowwise_instead_of_columnwise;
    
    public:

        explicit ConvergenceTable( std::string table_name = "---------- Default Table Name ----------" );
        
        ConvergenceTable& operator<<( Float entry );
        
        ConvergenceTable& operator<<( const std::string& seriesheader );
        
        
        ConvergenceTable& operator<<( char code );
        
        std::string text() const;
        
        std::string text( bool display_convergence_rates ) const;
        
        void lg() const;
        
        void lg( bool display_convergence_rates ) const;
        
        void print( std::ostream& os );
        


        void print( std::ostream& os, bool display_convergence_rates ) const;
        
        
        
        void print_standard( std::ostream& os, bool display_convergence_rates ) const;


        void print_transpose( std::ostream& os, bool display_convergence_rates ) const;

};




#endif
