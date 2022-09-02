
#include <cmath>

#include <iostream>
#include <chrono>

#include "../../basic.hpp"


using namespace std;

int main()
{
    cout << "Unit Test and Benchmark for Factorials and Binomials" << nl;
    
    /* survey factorials */

    std::cout << "\nLargest number whose factorial fits into data type (signed)" << std::endl;
    std::cout << "    case signed char       : " << largest_factorial_base<       signed char>() << std::endl;
    std::cout << "    case signed short      : " << largest_factorial_base<      signed short>() << std::endl;
    std::cout << "    case signed int        : " << largest_factorial_base<        signed int>() << std::endl;
    std::cout << "    case signed long       : " << largest_factorial_base<       signed long>() << std::endl;
    std::cout << "    case signed long long  : " << largest_factorial_base<  signed long long>() << std::endl;

    std::cout << "\nLargest number whose factorial fits into data type (unsigned)" << std::endl;
    std::cout << "    case unsigned char     : " << largest_factorial_base<     unsigned char>() << std::endl;
    std::cout << "    case unsigned short    : " << largest_factorial_base<    unsigned short>() << std::endl;
    std::cout << "    case unsigned int      : " << largest_factorial_base<      unsigned int>() << std::endl;
    std::cout << "    case unsigned long     : " << largest_factorial_base<     unsigned long>() << std::endl;
    std::cout << "    case unsigned long long: " << largest_factorial_base<unsigned long long>() << std::endl;

    /* check factorials */

    std::cout << "\nChecking factorial implementations...\n";

    for( int n = 0; n < 50; n++ ) {
        auto x = largest_factorial_base_AUX( n, 2 );
        // printf( "Largest k such that k! <= %d : %lld\n", n, (long long int)x );
        if(  0 <= n && n <=   1 ) Assert( x == 1 );
        if(  2 <= n && n <=   5 ) Assert( x == 2 );
        if(  6 <= n && n <=  23 ) Assert( x == 3 );
        if( 24 <= n && n <= 124 ) Assert( x == 4 );
    }

    
    for( int i = 0; i <= largest_factorial_base<int>(); i++ ) {

        int64_t f1 = factorial_integer( i );
        int64_t f2 = factorial_integer_naive( i );
        int64_t f3 = factorial_integer_table( i );
        int64_t f4 = factorial_integer_loop( i );

        assert( f1 == f2 and f2 == f3 and f2 == f4 );

    }

    for( int i = 0; i < 20; i++ ) {

        Float f1 = factorial_numerical( i );
        Float f2 = factorial_numerical_naive( i );
        Float f3 = factorial_numerical_table( i );
        Float f4 = factorial_numerical_loop( i );

        assert( isaboutequal( f1/f2, 1. ) );
        assert( isaboutequal( f1/f3, 1. ) );
        assert( isaboutequal( f1/f4, 1. ) );

        assert( isaboutequal( f1, f2 ) );
        assert( isaboutequal( f1, f3 ) );
        assert( isaboutequal( f1, f4 ) );

    }


    /* TODO: check binomials */


    /* Benchmark factorial computation */

    auto N = 400000000;
    auto M = 20;

    std::cout << "\nBenchmarking methods for computing " << N << " factorials...\n";

    srand(0);
    auto start_naive = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_integer_naive( rand() % M );
    auto end_naive   = std::chrono::system_clock::now();

    srand(0);
    auto start_loop = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_integer_loop( rand() % M );
    auto end_loop   = std::chrono::system_clock::now();

    srand(0);
    auto start_table = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_integer_table( rand() % M );
    auto end_table   = std::chrono::system_clock::now();

    srand(0);
    auto start_n_naive = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_numerical_naive( rand() % M );
    auto end_n_naive   = std::chrono::system_clock::now();

    srand(0);
    auto start_n_loop = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_numerical_loop( rand() % M  );
    auto end_n_loop   = std::chrono::system_clock::now();

    srand(0);
    auto start_n_table = std::chrono::system_clock::now();
    for( int i = 0; i < N; i++ ) factorial_numerical_table( rand() % M );
    auto end_n_table   = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds_naive   = end_naive   - start_naive;
    std::chrono::duration<double> elapsed_seconds_loop    = end_loop    - start_loop;
    std::chrono::duration<double> elapsed_seconds_table   = end_table   - start_table;
    std::chrono::duration<double> elapsed_seconds_n_naive = end_n_naive - start_n_naive;
    std::chrono::duration<double> elapsed_seconds_n_loop  = end_n_loop  - start_n_loop;
    std::chrono::duration<double> elapsed_seconds_n_table = end_n_table - start_n_table;
    
    std::cout << "Integer, Naive method, elapsed time:   " << elapsed_seconds_naive.count()   << "s\n";
    std::cout << "Integer, Loop method, elapsed time:    " << elapsed_seconds_loop.count()    << "s\n";
    std::cout << "Integer, Table method, elapsed time:   " << elapsed_seconds_table.count()   << "s\n";
    std::cout << "Numerical, Naive method, elapsed time: " << elapsed_seconds_n_naive.count() << "s\n";    
    std::cout << "Numerical, Loop method, elapsed time:  " << elapsed_seconds_n_loop.count()  << "s\n";
    std::cout << "Numerical, Table method, elapsed time: " << elapsed_seconds_n_table.count() << "s\n";    
    
    cout << "Finished Unit Test" << nl;
    
    return 0;
}
