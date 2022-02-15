
    #include <chrono>
    #include <cmath>
    #include <ctime>
    #include <cstdio>
    #include <vector>


    #define RUN \
            double result = 0.; \
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now(); \
            for( int l = 0; l < repetitions; l++ ) { \
                for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100; \
                for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i]; \
            } \
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now(); \
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count(); 


    int main(){

        int arraylength = 1000; // reduce this number in case your program crashes
        int repetitions = 200;
        
        printf( "array length = %d repetitions = %d\n", arraylength, repetitions );
        
        {
            srand( 1 + std::time(nullptr) );
            double data[arraylength]; // technically, non-compliant with C++ standard.
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < repetitions; l++ ) {
                for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            printf( "Calculation result: %e\n", result );
            const char * name = "C99 style stack array";
            printf( "Duration of %s: %ldms\n", name, duration );
        }

        {
            srand( 2 + std::time(nullptr) );
            std::vector<double> data( arraylength );
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < repetitions; l++ ) {
                for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            printf( "Calculation result: %e\n", result );
            const char * name = "std::vector array";
            printf( "Duration of %s: %ldms\n", name, duration );
        }

        {
            srand( std::time(nullptr) );
            double * data = new double[arraylength];
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < repetitions; l++ ) {
                for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            printf( "Calculation result: %e\n", result );
            const char * name = "C++ style heap array";
            printf( "Duration of %s: %ldms\n", name, duration );
            delete[] data;
        }

        {
            srand( std::time(nullptr) );
            double * data = (double*)malloc( arraylength * sizeof(double) );
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < repetitions; l++ ) {
                for( int i = 0; i < arraylength; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < arraylength; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            printf( "Calculation result: %e\n", result );
            const char * name = "C style heap array";
            printf( "Duration of %s: %ldms\n", name, duration );
            free(data);
        }

        return 0;
    }
