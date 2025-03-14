#include <iostream>
#include <cstdio>

#include "../../base/include.hpp"


int main() {

    // Initialization
    safedouble fromDouble(3.14);
    safedouble fromFloat(3.14f);
    safedouble fromLongDouble(3.14159265358979323846L);

    // Explicit conversion to double
    double d1 = static_cast<double>(fromDouble);
    double d2 = static_cast<double>(fromFloat);
    double d3 = static_cast<double>(fromLongDouble);

    // Using with printf (requires explicit cast)
    printf("safedouble as double: %.15f\n", static_cast<double>(fromDouble));
    printf("safedouble as double: %.15f\n", static_cast<double>(fromFloat));
    printf("safedouble as double: %.15f\n", static_cast<double>(fromLongDouble));

    // Using with printf (requires explicit cast)
    printf("safedouble as double, C-style cast: %.15f\n", (double)(fromDouble));
    printf("safedouble as double, C-style cast: %.15f\n", (double)(fromFloat));
    printf("safedouble as double, C-style cast: %.15f\n", (double)(fromLongDouble));

    // Using with iostreams
    std::cout << "safedouble with std::cout: " << d1 << "\n";
    std::cout << "safedouble with std::cout: " << d2 << "\n";
    std::cout << "safedouble with std::cout: " << d3 << "\n";

    // Example of an overflow
    // safedouble fromTooLargeLongDouble(1e310L); // This will throw

    return 0;
}
