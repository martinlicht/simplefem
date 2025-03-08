#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING
bool log_has_a_fresh_line = true;

#include <cmath>
#include <cstdio>
#include <ctime>

#include <limits> // for std::float_round_style

// For controlling the floating-point behavior
#include <cfenv>
#include <xmmintrin.h>
#include <cfloat>



enum class TextColors : unsigned char
{
    black          = 30,
    black_light    = 90,
    
    red            = 31,
    red_light      = 91,
    
    green          = 32,
    green_light    = 92, 
    
    yellow         = 33, 
    yellow_bright  = 93, 
    
    blue           = 34,
    blue_bright    = 94,
    
    magenta        = 35,
    magenta_bright = 95,
    
    cyan           = 36,
    cyan_light     = 96, 
    
    white          = 37,
    white_light    = 97 
};


Logger::Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline,
            const char* filename,
            const int linenumber
        )
        : 
        use_cerr( use_cerr ), // internalstream( use_cerr ? std::cerr : std::cout ), // internalstream( os ),
        pad_newline_if_there_is_none( do_newline ),
        filename( filename ),
        linenumber( linenumber )
        {
            // *this << std::setprecision(10);
        }


// auto digitalcodenow() -> std::string;

Logger::~Logger() noexcept
{
    
    FILE* f = ( use_cerr ? stderr : stdout );
    
    const auto& message_string = internal;
    
    
    // if the log string is empty, just do not print anything 

    if( message_string.empty() ) return;


    // complete the prefix string 

    const std::string time_code = digitalcodenow(); 

    #ifdef USE_COLORED_OUTPUT
    const std::string colorcode_begin = "\033[96m";
    const std::string colorcode_close = "\033[m";
    #else
    const std::string colorcode_begin = "";
    const std::string colorcode_close = "";
    #endif
    const std::string prefix = printf_into_string(
        "%s[%s %s %4d]%s\t", 
        colorcode_begin.c_str(), time_code.c_str(), filename.c_str(), linenumber, colorcode_close.c_str() );
    
    
    // assemble final output string 
    // reserve enough space so that (in principle) each nl can be prefixed

    int number_of_newlines = 0;
    for( const char character : message_string ) 
        if( character == '\n' )
            number_of_newlines++;

    std::string output_string;

    output_string.reserve( 1 + message_string.size() + number_of_newlines * prefix.size() );


    // if we have inherited a fresh line, then insert a prefix 
    
    if( log_has_a_fresh_line ){
        log_has_a_fresh_line = false;
        output_string.insert(0,prefix);
    }

    
    // for each nl (except possibly at the final position) insert a prefix after it
    
    for( int c = 0; c < message_string.length(); c++ ) {
    
        const char character = message_string[c];

        output_string += character;
        // fputc( character, f );

        if( character == '\n' and c != message_string.length()-1 ) 
            output_string += prefix;
        
    }
    
    
    // if the last character is not a newline and we automatically append a newline,
    // then modify the string accordingly
    if( message_string.back() != nl and pad_newline_if_there_is_none )
        output_string += nl;

    
    log_has_a_fresh_line = ( output_string.back() == nl );
    
    // if( message_string.back() == '\n' or pad_newline_if_there_is_none ) {
    //    
    //     log_has_a_fresh_line = true;
    //  
    //     if( message_string.back() != '\n' and pad_newline_if_there_is_none ) 
    //         output_string += nl;
    //         // fputc( nl, f ); 
    //
    // }

    fputs( output_string.c_str(), f );
    fflush(f);

}






#if defined(_OPENMP)
#include <omp.h>
#endif // #if defined(_OPENMP)

System_Reporter::System_Reporter() noexcept
{
    output();
}

System_Reporter::~System_Reporter() noexcept
{
    #if defined(_OPENMP)
    // LOG << "###\tSystem Reporter: finished\n";
    #endif
}   

void System_Reporter::output()
{
    
    #if defined(GIT_COMMIT_ID)
    LOG << "###\tCurrent Git commit ID: " << GIT_COMMIT_ID << nl;
    #endif

    const time_t t = std::time(nullptr);
    const tm zeit = *std::localtime(&t);
    LOGPRINTF("###\t%d-%02d-%02d %02d:%02d:%02d\n", zeit.tm_year + 1900, zeit.tm_mon + 1, zeit.tm_mday, zeit.tm_hour, zeit.tm_min, zeit.tm_sec );
    
    #if defined(_OPENMP)
    LOG << "###\tOpenMP Value: " << _OPENMP << nl;
    LOG << "###\tMaximum number of threads: " << omp_get_max_threads() << nl;
    LOG << "###\tMaximum active levels: " << omp_get_max_active_levels() << nl;
    // LOG << "###\tNested parallelism supported: " << ( omp_get_nested() ? "Yes" : "No" ) << nl;
    LOG << "###\tDynamic adjustment of the number of threads enabled: " << ( omp_get_dynamic() ? "Yes" : "No" ) << nl;
    LOG << "###\tThread limit: " << omp_get_thread_limit() << nl;
    LOG << "###\tNumber of processors available: " << omp_get_num_procs() << nl;
    LOG << "###\tDefault thread affinity: " << ( omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity" ) << nl;
    // LOG << "###\tMaximum number of processors: " << p << " -> " << omp_get_place_num_procs() << nl;
    LOG << "###\tMaximum number of places: " << omp_get_num_places() << nl;
    for( int p = 0; p < omp_get_num_places(); p++ ) {
        LOG << "###\t\tMaximum number of processors per place: " << p << " - > " << omp_get_place_num_procs(p) << nl;
    }
    #else
    LOG << "###\tOpenMP is not enabled.\n";
    #endif

    // #ifdef _OPENMP
    // LOGPRINTF("OpenMP is enabled.\n");
    // LOGPRINTF("Number of threads: %d\n", omp_get_max_threads());
    // LOGPRINTF("Nested parallelism supported: %s\n", omp_get_nested() ? "Yes" : "No");
    // LOGPRINTF("Dynamic adjustment of the number of threads enabled: %s\n", omp_get_dynamic() ? "Yes" : "No");
    // LOGPRINTF("Number of processors available: %d\n", omp_get_num_procs());
    // LOGPRINTF("Default thread affinity: %s\n", omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity");
    // #else
    // LOGPRINTF("OpenMP is not enabled.\n");
    // #endif

    LOGPRINTF("###\tCompiler version: %s\n", __VERSION__);

    #ifdef _DEBUG
    LOGPRINTF("###\tMSVC Debugging flags enabled.\n");
    #else
    LOGPRINTF("###\tMSVC Debugging flags not enabled.\n");
    #endif

    #ifdef NDEBUG
    LOGPRINTF("###\tNDEBUG enabled.\n");
    #else
    LOGPRINTF("###\tNDEBUG not enabled.\n");
    #endif

    #ifdef __OPTIMIZE__
    LOGPRINTF("###\tCompiler optimization level: %d\n", __OPTIMIZE__);
    // Add more performance-related information here
    #endif

    #ifdef __unix__
    LOGPRINTF("###\tOS: Unix\n");
    #elif defined(_WIN32)
    LOGPRINTF("###\tOS: Windows\n");
    #else
    #error "Unknown platform"
    #endif





    LOG << "###\tFLT_EVAL_METHOD: " << FLT_EVAL_METHOD << space;

    switch (FLT_EVAL_METHOD) {
        case -1:
            LOG << "The default precision is not known." << nl;
            break;
        case 0:
            LOG << "All operations and constants evaluate in the range and precision of the type used." << nl;
            break;
        case 1:
            LOG << "All operations and constants evaluate in the range and precision of double." << nl;
            break;
        case 2:
            LOG << "All operations and constants evaluate in the range and precision of long double." << nl;
            break;
        default:
            if (FLT_EVAL_METHOD < 0) {
                LOG << "Implementation-defined behavior." << nl;
            } else {
                LOG << "Unknown evaluation mode." << nl;
            }
            break;
    }


    LOG << "###\tFLT_ROUNDS: " << FLT_ROUNDS << space;

    switch (FLT_ROUNDS) {
        case std::round_indeterminate:
            LOG << "Rounding style cannot be determined." << nl;
            break;
        case std::round_toward_zero:
            LOG << "Rounding toward zero." << nl;
            break;
        case std::round_to_nearest:
            LOG << "Rounding toward nearest representable value." << nl;
            break;
        case std::round_toward_infinity:
            LOG << "Rounding toward positive infinity." << nl;
            break;
        case std::round_toward_neg_infinity:
            LOG << "Rounding toward negative infinity." << nl;
            break;
        default:
            LOG << "Unknown rounding style." << nl;
            break;
    }


    #if __cplusplus >= 201703L
    LOG << "###\tFLT_HAS_SUBNORM: " << FLT_HAS_SUBNORM << space;
    switch (FLT_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in float is indeterminable." << nl;
            break;
        case 0:
            LOG << "Float type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Float type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for FLT_HAS_SUBNORM." << nl;
            break;
    }

    LOG << "###\tDBL_HAS_SUBNORM: " << DBL_HAS_SUBNORM << space;
    switch (DBL_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in double is indeterminable." << nl;
            break;
        case 0:
            LOG << "Double type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Double type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for DBL_HAS_SUBNORM." << nl;
            break;
    }

    LOG << "###\tLDBL_HAS_SUBNORM: " << LDBL_HAS_SUBNORM << space;
    switch (LDBL_HAS_SUBNORM) {
        case -1:
            LOG << "Support for subnormal numbers in long double is indeterminable." << nl;
            break;
        case 0:
            LOG << "Long double type does not support subnormal numbers." << nl;
            break;
        case 1:
            LOG << "Long double type supports subnormal numbers." << nl;
            break;
        default:
            LOG << "Unknown value for LDBL_HAS_SUBNORM." << nl;
            break;
    }
    #endif 

    
    #if true and defined(_WIN32)
    {
        // Query the current floating-point control word
        unsigned int currentControlWord = _controlfp(0, 0);
        
        // Check precision control bits
        switch(currentControlWord & _MCW_PC) {
            case _PC_24: 
                LOG << "###\tPrecision control: 24-bit." << nl;
                break;
            case _PC_53:
                LOG << "###\tPrecision control: 53-bit." << nl;
                break;
            case _PC_64:
                LOG << "###\tPrecision control: 64-bit." << nl;
                break;
            default:
                LOG << "###\tPrecision control: unknown." << nl;
                break;
        }
    }    
    #endif

    #if true and defined(_WIN32)
    LOGPRINTF("###\tFlushing subnormal numbers\n");
    _controlfp_s( nullptr, _DN_FLUSH, _MCW_DN ); // Flush denormals, both operands and results, to zero 
    #elif false and defined(__SSE__) // TODO: fix this stuff in gcc
    LOGPRINTF("###\tFlushing subnormal numbers\n"); // ??? xmmintrin.h
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
    #endif

}

const System_Reporter system_reporter;




#endif
