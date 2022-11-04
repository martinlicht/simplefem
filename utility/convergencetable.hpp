#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

// #include <ostream>
#include <string>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    typedef long double Float; 
    
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
        
        void insert_numerical_entry( Float entry );
        void insert_seriesheader( const std::string& seriesheader );
        void insert_newline();


        std::string text() const;
        
        void lg() const;
        
        void lg( bool display_convergence_rates ) const;
        
//         void print( std::ostream& os );
        

        std::string text( bool display_convergence_rates ) const;
        
        std::string text_standard( bool display_convergence_rates ) const;

        std::string text_transpose( bool display_convergence_rates ) const;

};




ConvergenceTable& operator<<( ConvergenceTable& contable, float entry );

ConvergenceTable& operator<<( ConvergenceTable& contable, double entry );

ConvergenceTable& operator<<( ConvergenceTable& contable, long double entry );

ConvergenceTable& operator<<( ConvergenceTable& contable, const std::string& seriesheader );

ConvergenceTable& operator<<( ConvergenceTable& contable, char code );
        


#endif
