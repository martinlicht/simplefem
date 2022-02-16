
#include "../../basic.hpp"
#include "../../utility/utility.hpp"
#include <ostream>

using namespace std;

int main()
{
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
        
        return 0;
}
