
#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include <iostream>

using namespace std;

int main()
{
        ConvergenceTable Contable;
        
        Contable.lg();
        
        for( int i = 0; i < 5; i++ ) {
            for( int j = 0; j < 5; j++ )
                Contable << (Float)10. * i + j + 1;
            Contable << nl;
        }
        
        Contable.lg();
        
        return 0;
}
