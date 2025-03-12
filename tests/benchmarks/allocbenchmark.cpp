
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>


#define RUN \
        srand( seed ); \
        double result = 0.; \
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); \
        for( int l = 0; l < repetitions; l++ ) { \
            for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100; \
            for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i]; \
        } \
        std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now(); \
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count(); \
        printf( "Sum: %26.20e\t%10ju Microseconds\t", result, (uintmax_t)duration ); \
        printf( "%s\n", name )

        // printf( "Sum: %10.10e\t%10ld%lcs\t", result, duration, L'\u00b5' )


int main( int argc, char *argv[] )
{
    LOG << "Unit Test: allocation benchmark" << nl;
    
    int arraylength = 1000; // reduce this number in case your program crashes
    int repetitions = 800000;
    
    printf( "array length = %d\n", arraylength );
    printf( "repetitions  = %d\n", repetitions );
    
    auto seed = 1 + std::time(nullptr);

    {
        const char * name = "C99 style stack array";
        double data[arraylength]; // technically, non-compliant with C++ standard.
        RUN;
    }

    {
        const char * name = "std::vector array";
        std::vector<double> data( arraylength );
        RUN;
    }

    {
        const char * name = "C style heap array";
        double * data = (double*)malloc( arraylength * sizeof(double) );
        RUN;
        free(data);
    }

    {
        const char * name = "C++ style heap array";
        double * data = new double[arraylength];
        RUN;
        delete[] data;
    }

    return 0;
}
