#include <cstdio>
#include <cmath>
#include <cfloat>

#include <limits>

#include "../../basic.hpp"


template<typename NumericalClass>
struct SqrtTest
{
    
    
    struct sqrt_results{ 
        NumericalClass cpp; 
        NumericalClass fem; 
    };

    union numerical_container {
        NumericalClass f;
        unsigned int i;
    };

    sqrt_results calculate_sqrts( NumericalClass x ) {
        NumericalClass cpp = std::sqrt(x);
        NumericalClass fem = Sqrt(x,200);
        // std::printf( "%c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)x, (double)cpp, (double)fem );
        std::printf( "%c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)(safedouble)x, (double)(safedouble)cpp, (double)(safedouble)fem );
        // LOG << nl;
        // std::printf( "%10e %10e\n", (double)(cpp/fem), (double)(fem/cpp) );
        // LOG << ( std::isfinite(fem/cpp) ) << space << ( std::isfinite((double)(fem/cpp)) ) <<  nl;
        // LOG << ( std::isfinite(cpp/fem) ) << space << ( std::isfinite((double)(cpp/fem)) ) <<  nl;
        std::printf( "%10e %10e\n", (double)(safedouble)(cpp/fem), (double)(safedouble)(fem/cpp) );
        return { cpp, fem };
    }

    void check()
    {

        LOG << "\nSize of Numerical Type\n";

        {
            LOG << "SQRT: check machine epsilon" << nl;
            numerical_container u;
            u.f = std::numeric_limits<NumericalClass>::epsilon();
            calculate_sqrts( u.f );
        }

        { 
            LOG << "SQRT: very big integer" << nl;
            numerical_container u;
            u.f = std::numeric_limits<unsigned int>::max() >> 4;
            // LOG << u.i << ' ' << (double)u.f << nl;
            calculate_sqrts(u.f);
        }

        if(false)
        {
            LOG << "Test the ratio of sqrt implementations on range: " << nl;

            unsigned int iteration_min = 0;
            unsigned int iteration_max = 10000; 

            NumericalClass worst_upper = 0.;
            NumericalClass worst_lower = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(max:worst_upper,worst_lower)
            #endif
            for( int i = iteration_min; i < iteration_max; ++i )
            {
                numerical_container u;

                u.i = i;

                // Optionally, break or continue on specific conditions to avoid NaNs, infinities, etc.
                if (std::isnan(u.f) || std::isinf(u.f)) {
                    continue; // Skip this iteration
                }

                if( u.f < 0 ) continue; // skip from here on 

                auto sqrts = calculate_sqrts(u.f);

                auto ratio = sqrts.cpp / sqrts.fem;

                worst_upper = std::max( worst_upper, ratio );
                worst_lower = std::max( worst_lower, ratio );

                // LOGPRINTF( "%.10e %.10e %.10e \n", (double)(safedouble)worst_upper, (double)(safedouble)worst_lower );
            
            }

            LOGPRINTF( "WORST RATIOS: %e %e \n", (double)(safedouble)worst_upper, (double)(safedouble)worst_lower );
            
        }
            
        LOG << "Test finished\n";
        
    }


};




int main() 
{

    static_assert( sizeof(float) == sizeof(int), "sizeof(float) == sizeof(int) failed" );

    #if __cplusplus >= 201703L
    LOGPRINTF( "FLT_MIN=%e; TRUE_MIN=%e; HALF FLT_MIN=%e \n", FLT_MIN, FLT_TRUE_MIN, FLT_MIN / 2. );
    #endif

    SqrtTest<float> test_float;
    test_float.check();
    
    SqrtTest<double> test_double;
    test_double.check();
    
    SqrtTest<long double> test_longdouble;
    test_longdouble.check();
    
    

    return 0;
}
