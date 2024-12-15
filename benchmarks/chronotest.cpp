#include <iostream>
// #include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>

#if defined(__GNUC__) || defined(__clang__)
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

// the function f() does some time-consuming work
UNUSED static void f()
{
    volatile double d = 0;
    for(int n=0; n<10000; ++n)
       for(int m=0; m<10000; ++m)
           d = d + d*n*m;
}

volatile static long double foo;

int main( int argc, char *argv[] )
{
    
    std::clock_t c_start = std::clock();
    auto t_start = std::chrono::high_resolution_clock::now();
    
    // we implement the logistic map, just to kill some time. The value r is selected to ensure chaotic behavior
    // https://en.wikipedia.org/wiki/Logistic_map
    const long double r = 3.6;
    foo = 0.5;
    for( int t = 0; t < 1200000; t++) foo = r * foo * ( 1. - foo );
    
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
