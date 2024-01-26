#ifndef INCLUDEGUARD_STRINGS_HPP
#define INCLUDEGUARD_STRINGS_HPP

#ifdef FLAG_USE_CUSTOM_STRINGS

#include <string>

/////////////////////////////////////////////////
//                                             //
//      TEMPLATES TO TURN OBJECTS TO TEXT      //
//                                             //
/////////////////////////////////////////////////

inline std::string to_text( char input ) {
    return std::string(1,input);
}

inline std::string to_text( const char* input ) {
    return std::string(input);
}

inline std::string to_text( const std::string& input ) {
    return input;
}

inline std::string to_text( const std::string&& input ) {
    return input;
}

inline std::string to_text( const void* input ) {
    char buffer[ sizeof(decltype(input)) * 2 + 10 + 1 ]; // how pointers are printed is implementation-defined 
    std::snprintf( buffer, sizeof(buffer), "%p", input );
    return std::string(buffer);
}

template <typename T>
inline typename std::enable_if< std::is_integral<T>::value && std::is_signed<T>::value, std::string>::type
to_text(T input) {
    char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
    std::snprintf(buffer, sizeof(buffer), "%jd", static_cast<intmax_t>(input));
    return std::string(buffer);
}


template <typename T>
inline typename std::enable_if< std::is_integral<T>::value && std::is_unsigned<T>::value, std::string>::type
to_text(T input) {
    char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
    std::snprintf(buffer, sizeof(buffer), "%jd", static_cast<uintmax_t>(input));
    return std::string(buffer);
}

template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, std::string>::type
to_text( T input ) {
    char buffer[ std::numeric_limits<float>::max_digits10 + std::numeric_limits<float>::max_exponent10 + 10 + 1];
    std::snprintf( buffer, sizeof(buffer), "%Lg", (long double)input );
    return std::string(buffer);
}

// template <typename T, typename = decltype(std::declval<T>().text())>
// std::string to_text(const T& input) {
//     std::string text = input.text();
//     internal += text;
//     return *this;
// }

// std::string to_text( intmax_t input ) {
//     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
//     std::snprintf( buffer, sizeof(buffer), "%jd", input );
//     internal += buffer;
//     return *this;
// }

// std::string to_text( uintmax_t input ) {
//     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
//     std::snprintf( buffer, sizeof(buffer), "%ju", input );
//     internal += buffer;
//     return *this;
// }


#else // FLAG_USE_CUSTOM_STRINGS

#include <string>

template<typename T>
decltype( std::to_string( std::declval<T>() ) )
to_text( T t ) {
    return std::to_string( t );
}

#endif // FLAG_USE_CUSTOM_STRINGS




static_assert( std::is_same<int,signed int>::value, "int and signed int are not the same!" );
extern template std::string to_text( signed char );
extern template std::string to_text( signed short );
extern template std::string to_text( signed int );
extern template std::string to_text( signed long );
extern template std::string to_text( signed long long );
    
extern template std::string to_text( unsigned char );
extern template std::string to_text( unsigned short );
extern template std::string to_text( unsigned int );
extern template std::string to_text( unsigned long );
extern template std::string to_text( unsigned long long );

extern template std::string to_text( float );
extern template std::string to_text( double );
extern template std::string to_text( long double );





/////////////////////////////////////////////////
//                                             //
//            STRING UTILITIES                 //
//                                             //
/////////////////////////////////////////////////

/******************************************************/
/*      count the white space within STL string       */
/******************************************************/

int count_white_space( const std::string& str ); // TODO: Move to utilities 

/******************************************************/
/*          insert tabs before each line              */
/******************************************************/

std::string tab_each_line( std::string str ); // TODO: Move to utilities 









#endif
