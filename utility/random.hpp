#ifndef INCLUDEGUARD_UTILITY_RANDOM_HPP
#define INCLUDEGUARD_UTILITY_RANDOM_HPP


#include "../base/include.hpp"


// ============================================================================
// Utility functions for probability distributions
// ============================================================================

void seed_random_integer();

unsigned int random_integer();

unsigned int get_random_integer_modulo();

unsigned int flip_coin( Float prob_zero = 0.5 );

Float random_uniform();

Float gaussian_variable();



#endif
