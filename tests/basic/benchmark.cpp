
    #include <chrono>
    #include <cmath>
    #include <ctime>
    #include <iostream>
    #include <vector>


    int main(){

        int size = 1000; // reduce this number in case your program crashes
        int L = 200;
        
        std::cout << "size=" << size << " L=" << L << std::endl;
        {
            srand( 1 + std::time(nullptr) );
            double data[size]; // technically, non-compliant with C++ standard.
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < L; l++ ) {
                for( int i = 0; i < size; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < size; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "Calculation result is " << sqrt(result) << "\n";
            std::cout << "Duration of C99 style stack array: " << duration << "ms\n";
        }

        {
            srand( 2 + std::time(nullptr) );
            std::vector<double> data( size );
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < L; l++ ) {
                for( int i = 0; i < size; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < size; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "Calculation result is " << sqrt(result) << "\n";
            std::cout << "Duration of std::vector array:     " << duration << "ms\n";
        }

        {
            srand( std::time(nullptr) );
            double * data = new double[size];
            double result = 0.;
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for( int l = 0; l < L; l++ ) {
                for( int i = 0; i < size; i++ ) data[i] = rand() % 100;
                for( int i = 0; i < size; i++ ) result += data[i] * data[i];
            }
            std::chrono::steady_clock::time_point end   = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "Calculation result is " << sqrt(result) << "\n";
            std::cout << "Duration of C style heap array:    " << duration << "ms\n";
            delete[] data;
        }

        return 0;
    }
