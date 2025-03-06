
#include <cmath>
#include <limits>


#include "../basic/basic.hpp"
#include "random.hpp"

#include <random>


#ifdef _OPENMP 
static thread_local std::default_random_engine random_engine(0);
#else 
static              std::default_random_engine random_engine(0);
#endif
static thread_local std::uniform_int_distribution<unsigned int> distribution(0, std::numeric_limits<unsigned int>::max());

const constexpr unsigned int maximum_random_integer = std::numeric_limits<unsigned int>::max() - 10; // magic number

void seed_random_integer()
{
    // srand(0);
    random_engine.seed(0);
}

unsigned int random_integer()
{
    // int ret = rand();
    unsigned int ret = distribution( random_engine );
    
    Assert( 0 <= ret and ret <= std::numeric_limits<unsigned int>::max() );
    
    ret = ret % random_integer_maximum();
    
    return ret;
}

unsigned int random_integer_maximum()
{
    assert( maximum_random_integer > 0 );
    return maximum_random_integer;
}


unsigned int flip_coin( Float prob_zero )
{
    assert( 0 <= prob_zero and prob_zero <= 1. );
    Float value = random_uniform();
    if( value < prob_zero ) return 0; else return 1;
}

Float random_uniform()
{
    // Float ret = static_cast<Float>( rand() ) / static_cast<Float>( RAND_MAX );
    Float ret = static_cast<Float>( random_integer() ) / random_integer_maximum();
    Assert( 0. <= ret and ret <= 1. );
    return ret;
}


// Based on the implementations in the C-FAQ:
// http://c-faq.com/lib/gaussian.html

inline Float gaussrand_1()
{
    const int n = 25;
    
    Float x = 0;
    
    for( int i = 0; i < n; i++ ) 
        x += random_integer() / static_cast<Float>(random_integer_maximum());
    
    x -= n / 2.0;
    x /= std::sqrt( n / 12.0 );
    
    return x;
}

inline Float gaussrand_2()
{
    static bool phase = false;
    static Float u, v;
    const Float pi = 3.14159265358979323846;
    Float z;

    if( phase ) {
        z = std::sqrt( -2. * std::log(u) ) * std::cos( 2. * pi * v );
    } else {
        u = ( random_integer() + 1. ) / ( random_integer_maximum() + 2. );
        v = random_integer() / ( random_integer_maximum() + 1. );
        z = std::sqrt( -2. * std::log(u) ) * std::sin( 2. * pi * v );
    }
        
    phase = not phase;

    return z;
}

// http://c-faq.com/lib/gaussrand.luben.html
inline Float gaussrand_3( Float mean = 0., Float std_dev = 1. )
{
    Assert( std_dev > machine_epsilon );

    Float x = random_integer() / (random_integer_maximum() + 1.0);   /* 0.0 <= x < 1.0 */
    
    bool x_is_large = (x >= 0.5);
    
    Float y = std::abs(x - 1.0);                        /* 0.0 < y <= 1.0 */
    Float z = std_dev * std::sqrt( -2.0 * std::log(y) );

    return x_is_large ? (mean + z) : (mean - z);
}

inline Float gaussrand_4()
{
    const Float pi = 3.14159265358979323846;
    Float u = ( random_integer() + 1. ) / ( random_integer_maximum() + 2. );
    Float v = ( random_integer()      ) / ( random_integer_maximum() + 1. );

    Float result = std::sqrt( -2. * std::log(u) ) * std::sin( 2. * pi * v );

    return result;
}

Float gaussrand()
{
    Float ret = gaussrand_4();
    assert( std::isfinite(ret) );
    return ret;
}
