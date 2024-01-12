#include <iostream>
// #include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>
 
// the function f() does some time-consuming work
static void f()
{
    volatile double d = 0;
    for(int n=0; n<10000; ++n)
       for(int m=0; m<10000; ++m)
           d = d + d*n*m;
}

volatile long double foo;

int main()
{
    
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    foo = 1.;
    for( int t = 0; t < 1200000; t++) foo *= ( t - 1. ) / ( t + 1. );
    
    std::clock_t c_end = std::clock();
    auto t_end = std::chrono::high_resolution_clock::now();
 
    // std::cout << std::scientific << std::setprecision(2) << std::right;
    
    std::cout << "CPU time used: " << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
              << "CPU time used (unformatted): " << c_end-c_start << " ns\n"
              << "Wall clock time passed: " << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n";
              
    std::cout << "end = "
              // << std::setw(10)
              << c_end 
              << " and end = " << c_start << "ns\n"
              << "Wall clock time passed: " << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n";
    
}
