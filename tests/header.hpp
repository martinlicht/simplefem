#ifndef INCLUDEGUARD_UNITTEST_HEADER_HPP
#define INCLUDEGUARD_UNITTEST_HEADER_HPP

#include <strings>
#include <sstream>

#include "../basic.hpp"

ostream& log = std::cout;

#define TITLE    const title_string    = stringstream("") << 
#define SUBTITLE const subtitle_string = stringstream("") << 
#define NOTE     const note_string     = stringstream("") << 

#define MAIN int main() {

#define DONE return 0; }



#endif
