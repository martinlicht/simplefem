#ifndef INCLUDEGUARD_UTILITY_SUMMATION_HPP
#define INCLUDEGUARD_UTILITY_SUMMATION_HPP

#include <limits>

#include "../base/include.hpp"



template<typename NumericalType>
class NaiveSum {
private:
    NumericalType sum; // The running total

public:
    
    NaiveSum() : sum(0.0) {}
    void reset() { sum = 0.0; }
    NumericalType getSum() const { return sum; }
    operator NumericalType() const { return getSum(); }

    void add( NumericalType value ) { sum += value; }
    
    NumericalType operator+=( NumericalType value )
    {
        add( value ); 
        return getSum();
    }

};




template<typename NumericalType>
class KahanSum {
private:
    NumericalType sum;
    NumericalType compensation;

public:

    KahanSum() : sum(0.0), compensation(0.0) {}
    void reset() { sum = 0.0; compensation = 0.0; }
    NumericalType getSum() const { return sum; }
    operator NumericalType() const { return getSum(); }
    
    void add( NumericalType value ) {
        NumericalType y = value - compensation;// 1. Correct the value by removing the compensation
        NumericalType t = sum + y;             // 2. Add the corrected value to the running total
        compensation = (t - sum) - y;          // 3. Update the compensation for the next step
        sum = t;                               // 4. Update the running sum
    }

    NumericalType operator+=( NumericalType value )
    {
        add( value ); 
        return getSum();
    }
};





template<typename NumericalType>
class NeumaierSum {
private:
    NumericalType sum;
    NumericalType compensation;

public:
    
    NeumaierSum() : sum(0.0), compensation(0.0) {}
    void reset() { sum = 0.0; compensation = 0.0; }
    NumericalType getSum() const { return sum + compensation; }
    operator NumericalType() const { return getSum(); }
    
    void add( NumericalType value ) {
        NumericalType t = sum + value;         // 1. Add the value to the running total
        if( std::abs(sum) >= std::abs(value) )
            compensation += (sum - t) + value; // 2a. Update compensation if sum is larger
        else
            compensation += (value - t) + sum; // 2b. Update compensation if value is larger
        sum = t;                               // 3. Update the running total
    }

    NumericalType operator+=( NumericalType value )
    {
        add( value ); 
        return getSum();
    }
};




#endif // INCLUDEGUARD_UTILITY_SUMMATION_HPP
