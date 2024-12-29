#ifndef INCLUDEGUARD_UTILITY_SUMMATION_HPP
#define INCLUDEGUARD_UTILITY_SUMMATION_HPP

#include <limits>

#include "../basic.hpp"



template<typename NumericalType>
class NaiveSum {
private:
    NumericalType sum; // The running total

public:
    // Constructor to initialize the sum
    NaiveSum() : sum(0.0) {}

    // Add a value to the sum
    void add(NumericalType value) {
        sum += value; // Simply add the value to the running total
    }

    // Get the current total
    NumericalType getSum() const {
        return sum;
    }

    // Reset the summation
    void reset() {
        sum = 0.0;
    }

    NumericalType operator+=( NumericalType value )
    {
        add( value ); return getSum();
    }

    operator NumericalType() const 
    {
        return getSum();
    }
};




template<typename NumericalType>
class KahanSum {
private:
    NumericalType sum;          // The running total
    NumericalType compensation; // The compensation value

public:
    // Constructor to initialize the sum and compensation
    KahanSum() : sum(0.0), compensation(0.0) {}

    // Add a value to the sum
    void add(NumericalType value) {
        NumericalType y = value - compensation;       // Correct the value by removing the compensation
        NumericalType t = sum + y;                    // Add the corrected value to the running total
        compensation = (t - sum) - y;          // Update the compensation for the next step
        sum = t;                               // Update the running sum
    }

    // Get the current total
    NumericalType getSum() const {
        return sum;
    }

    // Reset the summation
    void reset() {
        sum = 0.0;
        compensation = 0.0;
    }

    NumericalType operator+=( NumericalType value )
    {
        add( value ); return getSum();
    }

    operator NumericalType() const 
    {
        return getSum();
    }
};





template<typename NumericalType>
class NeumaierSum {
private:
    NumericalType sum;          // The running total
    NumericalType compensation; // The compensation value

public:
    // Constructor to initialize the sum and compensation
    NeumaierSum() : sum(0.0), compensation(0.0) {}

    // Add a value to the sum
    void add(NumericalType value) {
        NumericalType t = sum + value;                // Add the value to the running total
        if (std::abs(sum) >= std::abs(value)) {
            compensation += (sum - t) + value; // Update compensation if sum is larger
        } else {
            compensation += (value - t) + sum; // Update compensation if value is larger
        }
        sum = t;                               // Update the running total
    }

    // Get the current total
    NumericalType getSum() const {
        return sum + compensation;            // Combine the running total and compensation
    }

    // Reset the summation
    void reset() {
        sum = 0.0;
        compensation = 0.0;
    }

    NumericalType operator+=( NumericalType value )
    {
        add( value ); return getSum();
    }

    operator NumericalType() const 
    {
        return getSum();
    }
};




#endif
