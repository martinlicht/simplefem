#ifndef INCLUDEGUARD_OPERATOR_IO_FLOATVECTOR_HPP
#define INCLUDEGUARD_OPERATOR_IO_FLOATVECTOR_HPP


#include <istream>
#include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"


/*******************
****  
****  Input / Output of float vectors 
****  
*******************/


void writeFloatVector( std::ostream& out, const FloatVector& vec );

FloatVector readFloatVector( std::istream& in );

void writeFloatVector( const char* filename, const FloatVector& vec );

FloatVector readFloatVector( const char* filename );



#endif
