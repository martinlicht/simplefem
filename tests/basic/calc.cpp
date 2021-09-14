
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "../../basic.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Factorials and Binomials" );

int main()
{
        LOG << "Unit Test: " << TestName << endl;
                
        cout << std::setprecision(10);

        assert( desired_precision < 1e-10 );

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





        auto N = 100000000;

        srand(0);
        auto start_naive = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_naive( rand() % 13 );
        auto end_naive   = std::chrono::system_clock::now();

        srand(0);
        auto start_loop = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_loop( rand() % 13  );
        auto end_loop   = std::chrono::system_clock::now();
    
        srand(0);
        auto start_table = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_table( rand() % 13  );
        auto end_table   = std::chrono::system_clock::now();

        srand(0);
        auto start_n_naive = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_numerical_naive( rand() % 13 );
        auto end_n_naive   = std::chrono::system_clock::now();

        srand(0);
        auto start_n_loop = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_numerical_loop( rand() % 13  );
        auto end_n_loop   = std::chrono::system_clock::now();


        std::chrono::duration<double> elapsed_seconds_naive = end_naive - start_naive;
        LOG << "\tNaive method, elapsed time:    " << elapsed_seconds_naive.count()   << "s\n";
        
        std::chrono::duration<double> elapsed_seconds_loop  = end_loop  - start_loop;
        LOG << "\tLoop method, elapsed time:     " << elapsed_seconds_loop.count()    << "s\n";
        
        std::chrono::duration<double> elapsed_seconds_table = end_table - start_table;
        LOG << "\tTable method, elapsed time:    " << elapsed_seconds_table.count()   << "s\n";
        
        std::chrono::duration<double> elapsed_seconds_n_naive = end_n_naive - start_n_naive;
        LOG << "\tNumerical, Naive method, elapsed time: " << elapsed_seconds_n_naive.count() << "s\n";
        
        std::chrono::duration<double> elapsed_seconds_n_loop  = end_n_loop  - start_n_loop;
        LOG << "\tNumerical, Loop method, elapsed time:  " << elapsed_seconds_n_loop.count()  << "s\n";
        
        
        LOG << "Finished Unit Test: " << TestName << endl;
        
        return 0;
}
