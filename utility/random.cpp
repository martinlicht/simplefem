
#include "../basic/basic.hpp"
#include "random.hpp"

void seed_random_integer()
{
    srand(0);
}

int random_integer()
{
    int ret = rand();
    Assert( 0 <= ret and ret <= RAND_MAX );
    return ret;
}

Float random_uniform()
{
    Float ret = static_cast<Float>( rand() ) / static_cast<Float>( RAND_MAX );
    Assert( 0. <= ret and ret <= 1. );
    return ret;
}


// Based on the implementations in the C-FAQ:
// http://c-faq.com/lib/gaussian.html

inline Float gaussrand_1()
{
    const int NSUM = 25;
    
    Float x = 0;
    
    for( int i = 0; i < NSUM; i++) 
        x += rand() / static_cast<Float>(RAND_MAX);
    
    x -= NSUM / 2.0;
    x /= std::sqrt( NSUM / 12.0 );
    
    return x;
}

inline Float gaussrand_2()
{
    static bool phase = false;
    static Float U, V;
    const Float PI = 3.14159265358979323846;
    Float Z;

    if( phase ) {
        Z = std::sqrt( -2. * std::log(U) ) * std::cos( 2. * PI * V );
    } else {
        U = ( rand() + 1. ) / ( RAND_MAX + 2. );
        V = rand() / ( RAND_MAX + 1. );
        Z = std::sqrt( -2. * std::log(U) ) * std::sin( 2. * PI * V );
    }
        
    phase = not phase;

    return Z;
}

// http://c-faq.com/lib/gaussrand.luben.html
inline Float gaussrand_3( Float mean = 0., Float std_dev = 1. )
{
    Assert( std_dev > machine_epsilon );

    Float x = rand() / (RAND_MAX + 1.0);   /* 0.0 <= x < 1.0 */
    
    bool large = (x < 0.5) ? false : true;
    
    Float y = std::abs(x - 1.0);                        /* 0.0 < y <= 1.0 */
    Float z = std_dev * std::sqrt( -2.0 * std::log(y) );

    return large ? (mean + z) : (mean - z);
}

inline Float gaussrand_4()
{
    const Float PI = 3.14159265358979323846;
    Float U = ( rand() + 1. ) / ( RAND_MAX + 2. );
    Float V = ( rand()      ) / ( RAND_MAX + 1. );
    return std::sqrt( -2. * std::log(U) ) * std::sin( 2. * PI * V );
}

Float gaussrand()
{
    return gaussrand_4();
}
