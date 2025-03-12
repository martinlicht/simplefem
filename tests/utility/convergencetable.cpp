
#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"

// using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test: convergence table" << nl;


    ConvergenceTable Contable("Test Table");
    
    Contable << "a";
    Contable << "b" << "c" << "d";
    Contable << "qwertyuiopasdfghjkl" << nl;

    Contable.lg();

    for( int i = 0; i < 7; i++ ) {
        for( int j = 0; j < 5; j++ )
            Contable << sign_power(j) * ( (Float)10. * i + j + 1 );
        Contable << nl;
    }

    // for( int i = 0; i < 5; i++ )
    //     Contable << -1.999 - i;
    // Contable << nl;
        
    // for( int i = 0; i < 7; i++ ) {
    //     for( int j = 0; j < 5; j++ )
    //         Contable << std::exp(-i*j);
    //     Contable << nl;
    // }

    assert( not Contable.print_rowwise_instead_of_columnwise );

    Contable.print_rowwise_instead_of_columnwise = false;        
    Contable.lg();
    
    Contable.print_rowwise_instead_of_columnwise = true;        
    Contable.lg();
    
    LOG << Contable.TeXtabular();
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
