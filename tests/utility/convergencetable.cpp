
#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include <iostream>

using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Convergence Table" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
    
        ConvergenceTable Contable;
        
        Contable.lg();
        
        for( int i = 0; i < 5; i++ ) {
            for( int j = 0; j < 5; j++ )
                Contable << (Float)10. * i + j + 1;
            Contable << nl;
        }
        
        Contable.lg();
        
        LOG << "Finished Unit Test: " << TestName << endl;
    
        return 0;
}
