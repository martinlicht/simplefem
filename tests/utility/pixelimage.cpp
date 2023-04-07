
#include "../../basic.hpp"
#include "../../utility/pixelimage.hpp"

using namespace std;

int main() {
    
    auto pim = readPixelImage("lena_color.tiff");

    for( size_t row = 0; row < pim.getheight(); row++ )
    for( size_t col = 0; col < pim.getwidth();  col++ )
    {
        pim(row,col).red   = 255 - pim(row,col).red;
        pim(row,col).green = 255 - pim(row,col).green;
        pim(row,col).blue  = 255 - pim(row,col).blue;
    }

    savePixelImage(pim,"anel_color.tiff");

    return 0;

}

