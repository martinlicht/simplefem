#ifndef INCLUDEGUARD_SPARSE_CSRRAINBOW_HPP
#define INCLUDEGUARD_SPARSE_CSRRAINBOW_HPP

//#include <memory>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "matcsr.hpp"



struct Rainbow
{
    Rainbow( const MatrixCSR& );

    void check() const;

    // number of colors used 
    int num_rows;
    int num_colors;

    // for each row index, give the color index 
    std::vector<int> F; 

    // CSR-like construction 
    std::vector<int> B; 
    std::vector<int> R;
};


#endif
