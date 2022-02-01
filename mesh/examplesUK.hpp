#ifndef INCLUDEGUARD_EXAMPLES_UNITEDKINGDOM_HPP
#define INCLUDEGUARD_EXAMPLES_UNITEDKINGDOM_HPP


#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "mesh.hpp"
#include "mesh.simplicial2D.hpp"


inline MeshSimplicial2D UnitedKingdom()
{
    // 191 triangles 
    
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 136, { // 136 points
        0,719,
        105,170,
        105,662,
        114,314,
        116,626,
        117,365,
        120,210,
        127,554,
        135,167,
        135,698,
        136,434,
        142,455,
        144,305,
        146,114,
        146,653,
        147,342,
        150,242,
        151,498,
        157,551,
        165,206,
        165,689,
        166,461,
        171,71,
        172,494,
        174,266,
        176,146,
        176,630,
        177,308,
        180,269,
        181,530,
        187,590,
        195,220,
        195,650,
        196,438,
        201,110,
        202,508,
        204,263,
        206,173,
        206,596,
        207,320,
        210,246,
        211,557,
        215,359,
        217,604,
        224,32,
        225,186,
        225,647,
        226,404,
        228,-24,
        231,124,
        232,474,
        234,302,
        236,150,
        236,608,
        237,359,
        240,212,
        241,534,
        24,288,
        245,398,
        247,570,
        254,71,
        255,162,
        255,686,
        256,416,
        258,14,
        261,90,
        262,450,
        264,316,
        266,116,
        266,647,
        267,358,
        270,224,
        271,500,
        27,324,
        275,412,
        277,546,
        284,70,
        285,188,
        285,700,
        286,455,
        288,28,
        291,66,
        292,476,
        294,282,
        296,128,
        296,646,
        297,319,
        300,263,
        301,512,
        305,378,
        307,572,
        314,31,
        315,220,
        315,666,
        316,454,
        321,92,
        322,508,
        -32,355,
        326,167,
        326,607,
        327,308,
        330,262,
        331,551,
        335,354,
        337,604,
        345,642,
        346,415,
        351,124,
        352,492,
        356,596,
        361,550,
        365,380,
        367,588,
        375,668,
        376,404,
        382,455,
        386,631,
        391,511,
        397,551,
        405,700,
        416,653,
        421,500,
        427,555,
        435,684,
        451,535,
        45,672,
        465,647,
        481,557,
        57,306,
        60,265,
        65,354,
        84,278,
        87,338,
        91,553,
        95,374,
        97,593,
      } ),
      {
        {100,86,103},
        {102,104,109},
        {102,109,110},
        {10,21,33},
        {102,90,104},
        {103,89,111},
        {104,105,109},
        {104,99,105},
        {105,113,116},
        {105,93,113},
        {106,94,115},
        {108,110,117},
        {108,96,110},
        {109,105,116},
        {110,109,112},
        {110,112,118},
        {111,106,114},
        {11,17,21},
        {112,109,116},
        {112,116,122},
        {113,119,120},
        {113,78,119},
        {114,106,115},
        {115,108,117},
        {116,113,120},
        {116,120,122},
        {117,110,118},
        {117,118,121},
        {118,112,122},
        {120,119,123},
        {121,118,122},
        {121,122,124},
        {12,15,27},
        {122,126,127},
        {124,122,127},
        {125,0,9},
        {128,130,132},
        {128,132,3},
        {128,73,130},
        {129,128,131},
        {131,12,16},
        {131,128,3},
        {132,130,134},
        {132,134,5},
        {13,25,34},
        {133,135,7},
        {135,4,7},
        {13,8,25},
        {14,20,32},
        {14,9,20},
        {1,6,8},
        {17,18,29},
        {18,30,41},
        {18,7,30},
        {19,31,37},
        {21,17,23},
        {21,23,50},
        {2,125,9},
        {22,34,44},
        {23,17,29},
        {23,29,35},
        {23,35,50},
        {24,12,27},
        {24,27,28},
        {24,28,40},
        {25,19,37},
        {25,37,52},
        {26,14,32},
        {28,27,36},
        {28,36,40},
        {2,9,14},
        {29,18,41},
        {30,26,38},
        {30,38,41},
        {31,24,40},
        {3,132,5},
        {31,40,55},
        {32,20,62},
        {33,50,66},
        {34,25,49},
        {34,49,60},
        {3,5,15},
        {35,29,41},
        {35,41,56},
        {36,27,39},
        {36,39,51},
        {37,31,45},
        {37,45,52},
        {38,26,43},
        {38,43,59},
        {40,36,51},
        {40,51,71},
        {41,38,59},
        {4,14,26},
        {4,2,14},
        {42,47,54},
        {43,32,46},
        {43,46,53},
        {43,53,59},
        {44,34,60},
        {44,60,64},
        {45,31,55},
        {45,55,61},
        {46,32,62},
        {46,62,69},
        {48,44,64},
        {49,25,52},
        {49,52,68},
        {50,35,56},
        {50,56,72},
        {51,39,54},
        {51,54,67},
        {51,67,83},
        {52,45,61},
        {52,61,68},
        {53,46,69},
        {53,69,90},
        {54,47,58},
        {54,58,70},
        {55,40,71},
        {55,71,77},
        {56,41,59},
        {56,59,75},
        {57,128,129},
        {57,73,128},
        {57,97,73},
        {58,63,74},
        {59,53,90},
        {60,76,80},
        {61,55,77},
        {6,16,19},
        {61,77,84},
        {62,20,78},
        {63,33,66},
        {63,66,74},
        {64,60,80},
        {65,49,68},
        {66,50,72},
        {66,72,82},
        {67,54,70},
        {67,70,86},
        {68,61,84},
        {69,62,78},
        {69,78,85},
        {69,85,90},
        {70,58,74},
        {70,74,89},
        {71,51,83},
        {71,83,87},
        {71,87,92},
        {72,56,75},
        {72,75,88},
        {7,4,30},
        {74,66,79},
        {74,79,94},
        {74,94,106},
        {75,59,90},
        {75,90,102},
        {76,68,84},
        {77,71,92},
        {79,66,82},
        {79,82,94},
        {80,76,81},
        {80,81,91},
        {8,19,25},
        {82,72,88},
        {82,88,96},
        {83,100,101},
        {83,67,86},
        {83,86,100},
        {84,77,98},
        {84,98,107},
        {85,78,93},
        {85,93,105},
        {86,70,89},
        {86,89,103},
        {87,83,101},
        {88,75,102},
        {89,106,111},
        {89,74,106},
        {9,0,78},
        {90,85,99},
        {90,99,104},
        {93,78,113},
        {94,108,115},
        {94,82,96},
        {94,96,108},
        {95,84,107},
        {96,102,110},
        {96,88,102},
        {99,85,105},
      }
    );
}





#endif
