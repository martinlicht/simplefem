#include <iostream>
#include <cassert>
#include <limits>
#include <cmath>

#include "../../basic.hpp"


int main( int argc, char *argv[] )
{

    LOG << "Test absolute\n";
    
    assert(absolute<int>(0) == 0);
    assert(absolute<int>(-10) == 10);
    assert(absolute<int>(10) == 10);
    assert(absolute<Float>(0.) == 0.);
    assert(absolute<Float>(-10.5) == 10.5);
    assert(absolute<Float>(10.5) == 10.5);
    assert(absolute<Float>( std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    assert(absolute<Float>(-std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    

    assert(absolute(0) == 0);
    assert(absolute(-5) == 5);
    assert(absolute(5) == 5);
    assert(absolute(0.) == 0.);
    assert(absolute(-3.14) == 3.14);
    assert(absolute(3.14) == 3.14);
    assert(absolute( std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    assert(absolute(-std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );

    LOG << "Test sign\n";
    
    assert(sign<int>(-10) == -1);
    assert(sign<int>(10) == 1);
    assert(sign<int>(0) == 0);
    assert(sign<Float>(-10.5) == -1);
    assert(sign<Float>(10.5) == 1);
    assert(sign<Float>(0.0) == 0);
    assert(sign<Float>( std::numeric_limits<Float>::infinity() ) ==  1 );
    assert(sign<Float>(-std::numeric_limits<Float>::infinity() ) == -1 );
    
    assert(sign(-5) == -1);
    assert(sign(5) == 1);
    assert(sign(0) == 0);
    assert(sign(-3.14) == -1.);
    assert(sign(3.14) == 1.);
    assert(sign(0.) == 0.);
    assert(sign( std::numeric_limits<Float>::infinity() ) ==  1 );
    assert(sign(-std::numeric_limits<Float>::infinity() ) == -1 );

    LOG << "Test sign_integer\n";

    assert(sign_integer<int>(-10) == -1);
    assert(sign_integer<int>(10) == 1);
    assert(sign_integer<int>(0) == 0);
    assert(sign_integer<Float>(-10.5) == -1);
    assert(sign_integer<Float>(10.5) == 1);
    assert(sign_integer<Float>(0.0) == 0);
    assert(sign_integer<Float>( std::numeric_limits<Float>::infinity() ) ==  1 );
    assert(sign_integer<Float>(-std::numeric_limits<Float>::infinity() ) == -1 );
    
    assert(sign_integer(-3) == -1);
    assert(sign_integer(3) == 1);
    assert(sign_integer(0) == 0);
    assert(sign_integer(-3.) == -1);
    assert(sign_integer(3.) == 1);
    assert(sign_integer(0.) == 0);
    assert(sign_integer( std::numeric_limits<Float>::infinity() ) ==  1 );
    assert(sign_integer(-std::numeric_limits<Float>::infinity() ) == -1 );

    LOG << "Test maximum\n";
    
    assert(maximum<int>(-7, -7) == -7);
    assert(maximum<int>( 0,  0) ==  0);
    assert(maximum<int>( 7,  7) ==  7);
    assert(maximum<Float>(-7.7, -7.7) == -7.7);
    assert(maximum<Float>( 0.0,  0.0) ==  0.0);
    assert(maximum<Float>( 7.7,  7.7) ==  7.7);
    
    assert(maximum<int>( 10,-10) == 10);
    assert(maximum<int>(  0,-10) ==  0);
    assert(maximum<int>(-10,-20) ==-10);
    assert(maximum<Float>(-10.5, 10.5) == 10.5);
    assert(maximum<Float>( 10.5,-10.5) == 10.5);
    assert(maximum<Float>(  0.5,-10.5) ==  0.5);
    assert(maximum<Float>(-10.5,-20.5) ==-10.5);
    assert(maximum<Float>( 1.0,  std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    assert(maximum<Float>( 1.0, -std::numeric_limits<Float>::infinity() ) == 1.0 );
    assert(maximum<Float>(  std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maximum<Float>( -std::numeric_limits<Float>::infinity(), 1.0 ) == 1.0 );
    
    assert(maximum(-10, 10) == 10);
    assert(maximum( 10,-10) == 10);
    assert(maximum(  0,-10) ==  0);
    assert(maximum( 10,-20) == 10);
    assert(maximum(-10.5, 10.5) == 10.5);
    assert(maximum( 10.5,-10.5) == 10.5);
    assert(maximum(  0.5,-10.5) ==  0.5);
    assert(maximum(-10.5,-20.5) ==-10.5);
    assert(maximum( 1.0, std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maximum( 1.0, -std::numeric_limits<Float>::infinity() ) == 1.0 );
    assert(maximum(  std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maximum( -std::numeric_limits<Float>::infinity(), 1.0 ) == 1.0 );
    assert(maximum(  std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maximum(  std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    
    LOG << "Test maxabs\n";
    
    assert(maxabs<int>(-7, -7) ==  7);
    assert(maxabs<int>( 0,  0) ==  0);
    assert(maxabs<int>( 7,  7) ==  7);
    assert(maxabs<Float>(-7.7, -7.7) ==  7.7);
    assert(maxabs<Float>( 0.0,  0.0) ==  0.0);
    assert(maxabs<Float>( 7.7,  7.7) ==  7.7);
    
    assert(maxabs<int>(-10, 10) == 10);
    assert(maxabs<int>( 10,-10) == 10);
    assert(maxabs<int>(  0,-10) == 10);
    assert(maxabs<int>(-10,-20) == 20);
    assert(maxabs<Float>(-10.5, 10.5) == 10.5);
    assert(maxabs<Float>( 10.5,-10.5) == 10.5);
    assert(maxabs<Float>(  0.5,-10.5) == 10.5);
    assert(maxabs<Float>(-10.5,-20.5) == 20.5);
    assert(maxabs<Float>( 1.0, std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maxabs<Float>( 1.0, -std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs<Float>(  std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs<Float>( -std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs<Float>( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maxabs<Float>( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    
    assert(maxabs(-10, 10) == 10);
    assert(maxabs( 10,-10) == 10);
    assert(maxabs(  0,-10) == 10);
    assert(maxabs(-10,-20) == 20);
    assert(maxabs(-10.5, 10.5) == 10.5);
    assert(maxabs( 10.5,-10.5) == 10.5);
    assert(maxabs(  0.5,-10.5) == 10.5);
    assert(maxabs(-10.5,-20.5) == 20.5);
    assert(maxabs( 1.0, std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maxabs( 1.0, -std::numeric_limits<Float>::infinity() ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs(  std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs( -std::numeric_limits<Float>::infinity(), 1.0 ) == std::numeric_limits<Float>::infinity() );
    assert(maxabs( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maxabs( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == std::numeric_limits<Float>::infinity() );
    assert(maxabs(-5, 10) == 10);
    assert(maxabs(-3.14, 2.71) == 3.14);

    LOG << "Test minimum\n";
    
    assert(minimum<int>(-7, -7) == -7);
    assert(minimum<int>( 0,  0) ==  0);
    assert(minimum<int>( 7,  7) ==  7);
    assert(minimum<Float>(-7.7, -7.7) == -7.7);
    assert(minimum<Float>( 0.0,  0.0) ==  0.0);
    assert(minimum<Float>( 7.7,  7.7) ==  7.7);
    
    assert(minimum<int>(-10, 10) ==-10);
    assert(minimum<int>( 10,-10) ==-10);
    assert(minimum<int>(  0,-10) ==-10);
    assert(minimum<int>(-10,-20) ==-20);
    assert(minimum<Float>(-10.5, 10.5) ==-10.5);
    assert(minimum<Float>( 10.5,-10.5) ==-10.5);
    assert(minimum<Float>(  0.5,-10.5) ==-10.5);
    assert(minimum<Float>(-10.5,-20.5) ==-20.5);
    assert(minimum<Float>( 1.0, std::numeric_limits<Float>::infinity() )  == 1.0 );
    assert(minimum<Float>( 1.0, -std::numeric_limits<Float>::infinity() ) == -std::numeric_limits<Float>::infinity() );
    assert(minimum<Float>(  std::numeric_limits<Float>::infinity(), 1.0 ) == 1.0 );
    assert(minimum<Float>( -std::numeric_limits<Float>::infinity(), 1.0 ) == -std::numeric_limits<Float>::infinity() );
    
    assert(minimum(-10, 10) ==-10);
    assert(minimum( 10,-10) ==-10);
    assert(minimum(  0,-10) ==-10);
    assert(minimum( 10,-20) ==-20);
    assert(minimum(-10.5, 10.5) ==-10.5);
    assert(minimum( 10.5,-10.5) ==-10.5);
    assert(minimum(  0.5,-10.5) ==-10.5);
    assert(minimum(-10.5,-20.5) ==-20.5);
    assert(minimum( 1.0, std::numeric_limits<Float>::infinity() )  == 1.0 );
    assert(minimum( 1.0, -std::numeric_limits<Float>::infinity() ) == -std::numeric_limits<Float>::infinity() );
    assert(minimum(  std::numeric_limits<Float>::infinity(), 1.0 ) == 1.0 );
    assert(minimum( -std::numeric_limits<Float>::infinity(), 1.0 ) == -std::numeric_limits<Float>::infinity() );
    assert(minimum( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == -std::numeric_limits<Float>::infinity() );
    assert(minimum( std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() )  == -std::numeric_limits<Float>::infinity() );
    assert(minimum(5, 10) == 5);
    assert(minimum(10, 5) == 5);
    assert(minimum(3.14, 2.71) == 2.71);

    
    
    
    LOG << "Test maximum, maxabs, minimum with multiple arguments or only one \n";

    LOG << "Test maximum with 1, 3, 4 arguments\n";
    
    assert(maximum<int>( 7 ) ==  7);
    assert(maximum<int>( 0 ) ==  0);
    assert(maximum<int>(-7 ) == -7);
    
    assert(maximum<int>( 7,  7,  7) ==  7);
    assert(maximum<int>( 0,  0,  0) ==  0);
    assert(maximum<int>(-7, -7, -7) == -7);
    
    assert(maximum<Float>( 7 ) ==  7);
    assert(maximum<Float>( 0 ) ==  0);
    assert(maximum<Float>(-7 ) == -7);
    
    assert(maximum<Float>( 7,  7,  7) ==  7);
    assert(maximum<Float>( 0,  0,  0) ==  0);
    assert(maximum<Float>(-7, -7, -7) == -7);

    assert(maximum<int>(10, 20, 5, 30) == 30);
    assert(maximum<Float>(-10.5, 10.5, 0.0) == 10.5);
    assert(maximum<Float>(-10.5, 30.0, 10.5, 0.0) == 30.0);

    LOG << "Test maxabs with 3, 4 arguments\n";

    assert(maxabs<int>( 7 ) ==  7);
    assert(maxabs<int>( 0 ) ==  0);
    assert(maxabs<int>(-7 ) ==  7);
    
    assert(maxabs<int>( 7,  7,  7) ==  7);
    assert(maxabs<int>( 0,  0,  0) ==  0);
    assert(maxabs<int>(-7, -7, -7) ==  7);
    
    assert(maxabs<Float>( 7 ) ==  7);
    assert(maxabs<Float>( 0 ) ==  0);
    assert(maxabs<Float>(-7 ) ==  7);
    
    assert(maxabs<Float>( 7,  7,  7) ==  7);
    assert(maxabs<Float>( 0,  0,  0) ==  0);
    assert(maxabs<Float>(-7, -7, -7) ==  7);

    assert(maxabs<int>(-10, 20, -5) == 20);
    assert(maxabs<int>(-10, 20, -5, 30) == 30);
    assert(maxabs<Float>(-10.5,  10.5, -5.0) == 10.5);
    assert(maxabs<Float>(-10.5,   4.5,  0.0) == 10.5);
    assert(maxabs<Float>(-10.5, -30.0,  1.5, 0.0) == 30.);

    LOG << "Test minimum with 3, 4 arguments\n";
    
    assert(minimum<int>( 7 ) ==  7);
    assert(minimum<int>( 0 ) ==  0);
    assert(minimum<int>(-7 ) == -7);
    
    assert(minimum<int>( 7,  7,  7) ==  7);
    assert(minimum<int>( 0,  0,  0) ==  0);
    assert(minimum<int>(-7, -7, -7) == -7);
    
    assert(minimum<Float>( 7 ) ==  7);
    assert(minimum<Float>( 0 ) ==  0);
    assert(minimum<Float>(-7 ) == -7);
    
    assert(minimum<Float>( 7,  7,  7) ==  7);
    assert(minimum<Float>( 0,  0,  0) ==  0);
    assert(minimum<Float>(-7, -7, -7) == -7);
    
    assert(minimum<int>(10, 20, 5) == 5);
    assert(minimum<int>(10, 20, 5, -30) == -30);
    assert(minimum<Float>(-10.5,  4.5, 0.0) == -10.5);
    assert(minimum<Float>( 10.5,-10.5, 0.0) == -10.5);
    assert(minimum<Float>(-10.5, -12.0, -10.0) == -12.0);




    LOG << "Test square\n";
    
    assert(square<int>(  0) == 0);
    assert(square<int>(-10) == 100);
    assert(square<int>( 10) == 100);
    assert(square<Float>(  0.) == 0.);
    assert(square<Float>(-10.5) == 10.5 * 10.5);
    assert(square<Float>(10.5) == 10.5 * 10.5);

    assert(square(0) == 0);
    assert(square(0.) == 0.);
    assert(square(5) == 25);
    assert(square(-5) == 25);
    assert(square(3.14) == 3.14 * 3.14);

    LOG << "Test is_numerically_small\n";
    
    assert( is_numerically_small(0.0, desired_closeness));
    assert(!is_numerically_small(1.0, desired_closeness));

    assert(not is_numerically_small( 2.0*desired_closeness ));
    assert(    is_numerically_small( 0.5*desired_closeness ));

    LOG << "Test is_numerically_close\n";
    
    assert( is_numerically_close(1.0, 1.0, desired_closeness));
    assert(!is_numerically_close(1.0, 2.0, desired_closeness));
    assert( is_numerically_close(  std::numeric_limits<Float>::infinity(),  std::numeric_limits<Float>::infinity() ) );
    assert( is_numerically_close( -std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() ) );
    Assert( not is_numerically_close(  std::numeric_limits<Float>::infinity(), -std::numeric_limits<Float>::infinity() ) );
    Assert( not is_numerically_close( -std::numeric_limits<Float>::infinity(),  std::numeric_limits<Float>::infinity() ) );

    assert( is_numerically_close( 1.0, 1.0 + desired_closeness/2. ));
    assert(!is_numerically_close(1.0, 1.1));

    LOG << "Test is_numerically_one\n";
    
    assert( is_numerically_one(1.0, desired_closeness));
    assert(!is_numerically_one(2.0, desired_closeness));

    assert( is_numerically_one(1.0));
    assert( is_numerically_one(1.0 + desired_closeness/2));
    assert(!is_numerically_one(1.1));



    {
        LOG << "Test power_numerical\n";
        assert(power_numerical(2.0,  0) == 1.0);
        assert(power_numerical(2.0,  1) == 2.0);
        assert(power_numerical(2.0,  3) == 8.0);
        assert(power_numerical(2.0, -2) == 0.25); 

        LOG << "Test power_integer\n";
        assert(power_integer( 2, 0) ==    1);
        assert(power_integer( 2, 1) ==    2);
        assert(power_integer( 2, 3) ==    8);
        assert(power_integer(-5, 0) ==    1);
        assert(power_integer(-5, 1) ==   -5);
        assert(power_integer(-5, 3) == -125);
        // assert(power_integer(0, 0) == 1); LOG << "Handle special case (0^0)

        LOG << "Test power_of_two\n";
        assert(power_of_two( 0) == 1);
        assert(power_of_two( 1) == 2);
        assert(power_of_two( 3) == 8);
        assert(power_of_two(10) == 1024);

        LOG << "Test sign_power\n";
        assert(sign_power(0) ==  1);
        assert(sign_power(1) == -1);
        assert(sign_power(2) ==  1);
        assert(sign_power(3) == -1);
    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    return 0;
}


