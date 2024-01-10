
#include "../../basic.hpp"
#include "../../utility/random.hpp"
#include "../../utility/pixelimage.hpp"

using namespace std;

int main( int argc, char** argv ) {
    
    

    // auto pim = readPixelImage("lena_color.tiff");
    // auto pim = readPixelImage("aurora.jpeg");
    // auto pim = readPixelImage("sanfrancisco.jpg");

    PixelImage pim( 300, 400 );

    LOG << pim.getheight() << "x" << pim.getwidth() << nl;

    for( size_t row = 0; row < pim.getheight(); row++ )
    for( size_t col = 0; col < pim.getwidth();  col++ )
    {
        pim(row,col).red   = random_integer() % 256;
        pim(row,col).green = random_integer() % 256;
        pim(row,col).blue  = random_integer() % 256;
    }

    savePixelImage(pim,"temp.tiff");

    for( size_t row = 0; row < pim.getheight(); row++ )
    for( size_t col = 0; col < pim.getwidth();  col++ )
    {
        pim(row,col).red   = 255 - pim(row,col).red;
        pim(row,col).green = 255 - pim(row,col).green;
        pim(row,col).blue  = 255 - pim(row,col).blue;
    }

    savePixelImage(pim,"reverse.tiff");

    return 0;

}

