#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

// #include <ostream>
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




ConvergenceTable& operator<<( ConvergenceTable& contable, float entry )
{
    contable.insert_numerical_entry( entry );
    return contable;
}

ConvergenceTable& operator<<( ConvergenceTable& contable, double entry )
{
    contable.insert_numerical_entry( entry );
    return contable;
}

ConvergenceTable& operator<<( ConvergenceTable& contable, long double entry )
{
    contable.insert_numerical_entry( entry );
    return contable;
}

ConvergenceTable& operator<<( ConvergenceTable& contable, const std::string& seriesheader )
{
    contable.insert_seriesheader( seriesheader );
    return contable;
}

ConvergenceTable& operator<<( ConvergenceTable& contable, char code )
{
    assert( code == nl );
    contable.insert_newline();
    return contable;
};
        


#endif
