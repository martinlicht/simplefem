#ifndef INCLUDEDGARD_BASIC_HPP
#define INCLUDEDGARD_BASIC_HPP

typedef double Float;

static const char space = ' ';

static const char nl = '\n';

static const char tab = '\t';

inline int kronecker( int i, int j )
{
    if( i == j )
        return 1;
    else
        return 0;
}



#endif
