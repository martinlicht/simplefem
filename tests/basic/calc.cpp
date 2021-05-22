
#include <cmath>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>

#include "../../basic.hpp"


using namespace std;

int main()
{
        cout << "Unit Test for Factorials and Binomials";// << endl;
        
        cout << std::setprecision(10);

        for( int i = 0; i <= 12; i++ ) {
            int64_t f1 = factorial_integer( i );
            int64_t f2 = factorial_integer_naive( i );
            int64_t f3 = factorial_integer_table( i );
            int64_t f4 = factorial_integer_loop( i );
            assert( f1 == f2 and f2 == f3 and f2 == f4 );
        }

        auto N = 100000000;

        srand(0);
        auto start_naive = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_naive( rand() % 13 );
        auto end_naive   = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds_naive = end_naive - start_naive;
        std::cout << "Naive method, elapsed time:    " << elapsed_seconds_naive.count()   << "s\n";
        
        srand(0);
        auto start_loop = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_loop( rand() % 13  );
        auto end_loop   = std::chrono::system_clock::now();
    
        std::chrono::duration<double> elapsed_seconds_loop  = end_loop  - start_loop;
        std::cout << "Loop method, elapsed time:     " << elapsed_seconds_loop.count()    << "s\n";
        
        srand(0);
        auto start_table = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_integer_table( rand() % 13  );
        auto end_table   = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds_table = end_table - start_table;
        std::cout << "Table method, elapsed time:    " << elapsed_seconds_table.count()   << "s\n";
        
        srand(0);
        auto start_n_naive = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_numerical_naive( rand() % 10 );
        auto end_n_naive   = std::chrono::system_clock::now();

        std::chrono::duration<double> elapsed_seconds_n_naive = end_n_naive - start_n_naive;
        std::cout << "N, Naive method, elapsed time: " << elapsed_seconds_n_naive.count() << "s\n";
        
        srand(0);
        auto start_n_loop = std::chrono::system_clock::now();
        for( int i = 0; i < N; i++ ) factorial_numerical_loop( rand() % 10  );
        auto end_n_loop   = std::chrono::system_clock::now();
    
        std::chrono::duration<double> elapsed_seconds_n_loop  = end_n_loop  - start_n_loop;
        std::cout << "N, Loop method, elapsed time:  " << elapsed_seconds_n_loop.count()  << "s\n";
        

        
        cout << "Finished Unit Test";// << endl;
        
        return 0;
}
