
#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include <ostream>

using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Convergence Table" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
    
        ConvergenceTable Contable("Test Table");
        
        Contable.lg();

        Contable << "a" << "b" << "c" << "d" << "qwertyuiopasdfghjkl";

        for( int i = 0; i < 7; i++ ) {
            for( int j = 0; j < 5; j++ )
                Contable << ( (Float)10. * i + j + 1 );
            Contable << nl;
        }

        for( int i = 0; i < 5; i++ )
            Contable << -1.999 - i;
        Contable << nl;
            
        for( int i = 0; i < 7; i++ ) {
            for( int j = 0; j < 5; j++ )
                Contable << std::exp(-i*j);
            Contable << nl;
        }

        Contable.lg();
        
        Contable.print_rowwise_instead_of_columnwise = true;
        
        Contable.lg();
        
        LOG << "Finished Unit Test: " << TestName << endl;
    
        return 0;
}
