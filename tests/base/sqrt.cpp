#include <cstdio>
#include <cmath>
#include <cfloat>

#include <limits>

#include "../../basic.hpp"


template<typename NumericalType>
struct SqrtTest
{
    struct SqrtResults{ 
        NumericalType cpp; 
        NumericalType fem; 
    };

    union numerical_container {
        NumericalType f;
        unsigned int i;
    };

    SqrtResults calculate_sqrts( NumericalType x ) {
        NumericalType cpp = std::sqrt(x);
        NumericalType fem = Sqrt(x);
        // LOGPRINTF( "result %c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)x, (double)cpp, (double)fem );
        // LOGPRINTF( "result %c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)(safedouble)x, (double)(safedouble)cpp, (double)(safedouble)fem );
        // LOG << nl;
        // LOGPRINTF( "%10e %10e\n", (double)(cpp/fem), (double)(fem/cpp) );
        // LOG << ( std::isfinite(fem/cpp) ) << space << ( std::isfinite((double)(fem/cpp)) ) <<  nl;
        // LOG << ( std::isfinite(cpp/fem) ) << space << ( std::isfinite((double)(cpp/fem)) ) <<  nl;
        // LOGPRINTF( "Ratios: %10e %10e\n", (double)(safedouble)(cpp/fem), (double)(safedouble)(fem/cpp) );
        return { cpp, fem };
    }

    void check()
    {

        // LOG << "\nSize of Numerical Type\n";
        LOG << nl;

        {
            LOG << "SQRT: check machine epsilon" << nl;
            numerical_container u;
            u.f = std::numeric_limits<NumericalType>::epsilon();
            auto result = calculate_sqrts( u.f );
            LOGPRINTF( "epsilon x=%.30e cpp=%.10e fem=%.10e \n", 
                (double)(safedouble)u.f, (double)(safedouble)result.cpp, (double)(safedouble)result.fem ); 
        }

        { 
            LOG << "SQRT: very big integer" << nl;
            numerical_container u;
            u.f = std::numeric_limits<unsigned int>::max() >> 6;
            // LOG << u.i << ' ' << (double)u.f << nl;
            auto result = calculate_sqrts( u.f );
            LOGPRINTF( "epsilon x=%.30e cpp=%.10e fem=%.10e \n", 
                (double)(safedouble)u.f, (double)(safedouble)result.cpp, (double)(safedouble)result.fem ); 
        }

        // if(false)
        {
            LOG << "Test the ratio of sqrt implementations on range: " << nl;

            unsigned int iteration_min = 0;
            unsigned int iteration_max = 10000; 

            NumericalType worst_upper = 0.;
            NumericalType worst_lower = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(max:worst_upper,worst_lower)
            #endif
            for( int i = iteration_min; i < iteration_max; ++i )
            {
                numerical_container u;

                u.f = i;

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

            LOGPRINTF( "WORST RATIOS: %.17e %.17e \n", (double)(safedouble)worst_upper, (double)(safedouble)worst_lower );
            
        }
            
        LOG << "Test finished\n\n";
        
    }


};















/*

template<typename NumericalType>
struct ThirdRootTest
{
    struct thirdroot_results{ 
        NumericalType cpp; 
        NumericalType fem; 
    };

    union numerical_container {
        NumericalType f;
        unsigned int i;
    };

    thirdroot_results calculate_thirdroots( NumericalType x ) {
        NumericalType cpp = std::pow( x, 1./3. );
        NumericalType fem = ThirdRoot(x);
        // LOGPRINTF( "result %c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)x, (double)cpp, (double)fem );
        // LOGPRINTF( "result %c x=%.30e cpp=%.10e fem=%.10e \n", std::isnormal(x) ? 'n' : 's', (double)(safedouble)x, (double)(safedouble)cpp, (double)(safedouble)fem );
        // LOG << nl;
        // LOGPRINTF( "%10e %10e\n", (double)(cpp/fem), (double)(fem/cpp) );
        // LOG << ( std::isfinite(fem/cpp) ) << space << ( std::isfinite((double)(fem/cpp)) ) <<  nl;
        // LOG << ( std::isfinite(cpp/fem) ) << space << ( std::isfinite((double)(cpp/fem)) ) <<  nl;
        // LOGPRINTF( "Ratios: %10e %10e\n", (double)(safedouble)(cpp/fem), (double)(safedouble)(fem/cpp) );
        return { cpp, fem };
    }

    void check()
    {

        // LOG << "\nSize of Numerical Type\n";

        LOG << nl;

        {
            LOG << "3RD: check machine epsilon" << nl;
            numerical_container u;
            u.f = std::numeric_limits<NumericalType>::epsilon();
            auto result = calculate_thirdroots( u.f );
            LOGPRINTF( "epsilon x=%.30e cpp=%.10e fem=%.10e \n", 
                (double)(safedouble)u.f, (double)(safedouble)result.cpp, (double)(safedouble)result.fem );         
        }

        { 
            LOG << "3RD: very big integer" << nl;
            numerical_container u;
            u.f = std::numeric_limits<unsigned int>::max() >> 6;
            // LOG << u.i << ' ' << (double)u.f << nl;
            auto result = calculate_thirdroots( u.f );
            LOGPRINTF( "epsilon x=%.30e cpp=%.10e fem=%.10e \n", 
                (double)(safedouble)u.f, (double)(safedouble)result.cpp, (double)(safedouble)result.fem ); 
        }

        // if(false)
        {
            LOG << "Test the ratio of third root implementations on range: " << nl;

            unsigned int iteration_min = 0;
            unsigned int iteration_max = 10000; 

            NumericalType worst_upper = 0.;
            NumericalType worst_lower = 0.;
            
            #if defined(_OPENMP)
            #pragma omp parallel for reduction(max:worst_upper,worst_lower)
            #endif
            for( int i = iteration_min; i < iteration_max; ++i )
            {
                numerical_container u;

                u.f = i;

                // Optionally, break or continue on specific conditions to avoid NaNs, infinities, etc.
                if (std::isnan(u.f) || std::isinf(u.f)) {
                    continue; // Skip this iteration
                }

                if( u.f < 0 ) continue; // skip from here on 

                auto thirdroots = calculate_thirdroots(u.f);

                auto ratio = thirdroots.cpp / thirdroots.fem;

                worst_upper = std::max( worst_upper, ratio );
                worst_lower = std::max( worst_lower, ratio );

                // LOGPRINTF( "%.10e %.10e %.10e \n", (double)(safedouble)worst_upper, (double)(safedouble)worst_lower );
            
            }

            LOGPRINTF( "WORST RATIOS: %.17e %.17e \n", (double)(safedouble)worst_upper, (double)(safedouble)worst_lower );
            
        }
            
        LOG << "Test finished\n\n";
        
    }


};

*/



int main() 
{

    static_assert( sizeof(float) == sizeof(int), "sizeof(float) == sizeof(int) failed" );

    #if __cplusplus >= 201703L
    LOGPRINTF( "FLT_MIN=%e; TRUE_MIN=%e; HALF FLT_MIN=%e \n", 
        std::numeric_limits<float>::min(),        // FLT_MIN 
        std::numeric_limits<float>::denorm_min(), // FLT_TRUE_MIN, 
        std::numeric_limits<float>::min() / 2.    // FLT_MIN / 2. 
    );
    #endif
    
    SqrtTest<float> test_sqrt_float;
    test_sqrt_float.check();
    
    SqrtTest<double> test_sqrt_double;
    test_sqrt_double.check();
    
    SqrtTest<long double> test_sqrt_longdouble;
    test_sqrt_longdouble.check();
    
    
/*
    ThirdRootTest<float> test_thirdroot_float;
    test_thirdroot_float.check();
    
    ThirdRootTest<double> test_thirdroot_double;
    test_thirdroot_double.check();
    
    ThirdRootTest<long double> test_thirdroot_longdouble;
    test_thirdroot_longdouble.check();
*/    
    

    return 0;
}
