#include <string>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef std::string RGBstring;

RGBstring rgb_to_string( unsigned char r, unsigned char g, unsigned char b )
{
    char result[8] = {0,0,0,0,0,0,0,0};
    sprintf(result, "#%02x%02x%02x", r, g, b);
    return std::string( result );
}

// unsigned char[3] string_to_rgb( RGBstring rgb )
// {
//     unsigned char r, g, b;
//     r = ( rgb[1] - 'A' ) * 16 + ( rgb[2] - 'A' );
//     g = ( rgb[3] - 'A' ) * 16 + ( rgb[4] - 'A' );
//     b = ( rgb[5] - 'A' ) * 16 + ( rgb[6] - 'A' );
//     return { r, g, b };
// }

int leading_digits( double num )
{
    // If the number is negative, take its absolute value
    // Calculate the order of magnitude of the number
    // If the number is greater than 1, return the integer part of the order, else, return 0
    
    if( num < 0 ) num = -num;
    int order = (int)ceil( log10(num) );
    if( order >= 0 ) return order + 1; else return 0;
}

std::string render_number( double num, int lead, int tail )
{
    char str[1+lead+1+tail+1];
    sprintf( str, "%0*.*f", lead+1+tail, tail, num);
    return std::string(str);
}







#include <Magick++.h>
#include <iostream>

using namespace Magick;

// g++ images.cpp `Magick++-config --cxxflags --cppflags --ldflags --libs`
// sudo apt install graphicsmagick-libmagick-dev-compat libgraphicsmagick++1-dev libmagick++-6-headers libmagick++-dev
// sudo apt-get update; sudo apt-get upgrade; sudo apt-get install imagemagick-common




/*
    <polygon points="10,20 40,50 0010,10" fill="#FF3399" stroke="none"></polygon>
    <polygon points="10,20 40,50 0010,10" fill="#FF3399" stroke="#11CCDD" stroke-width="1"></polygon>

    void writeSVG( 
        const std::string filename,
        const std::vector<std::array<3,int>>& triangles, 
        const Coordinates& coords,
        const FloatVector* triangle_values = nullptr
    ) {



    };

    int main() {
        unsigned char r = 255;
        unsigned char g = 0;
        unsigned char b = 128;
        RGBstring rgb = rgb_string( r, g, b );
        printf( "RGB string: %s\n", rgb );
        return 0;
    }

*/

enum class ColorComp { red, green, blue };

struct Pixel
{
    unsigned char red;
    unsigned char green;
    unsigned char blue;

    unsigned char operator()(ColorComp cc) const { 
        switch(cc) {
            case ColorComp::red:   return red;
            case ColorComp::green: return green;
            case ColorComp::blue:  return blue;
        }
    }
};

class PixelImage
{
    private:
        size_t height;
        size_t width;
        std::vector<Pixel> data;

    public:

        PixelImage( size_t height, size_t width ) 
        : height(height), width(width), data( width*height, {128,128,128} ) 
        {}
        
        ~PixelImage() {}

        size_t getheight(){ return height; };

        size_t getwidth(){ return width; };

        Pixel& operator()( size_t row, size_t col ) {
            assert(row >= 0 && row < height);
            assert(col >= 0 && col < width);
            return data.at( row * width + col );
        }
        
        Pixel operator()( size_t row, size_t col ) const {
            assert(row >= 0 && row < height);
            assert(col >= 0 && col < width);
            return data.at( row * width + col );
        }
};


PixelImage readPixelImage( std::string str )
{
    
    InitializeMagick(nullptr);

    Image image;    
    image.read( str.c_str() );

    size_t width  = image.columns();
    size_t height = image.rows();

    PixelImage ret( height, width ); 
    
    // Get the pixel values and store them
    Pixels view(image);
    const PixelPacket* pixel_cache = view.getConst(0, 0, width, height);
    for (size_t row = 0; row < height; row++) 
    for (size_t col = 0; col < width; col++) 
    {
            size_t index = row * width + col;

            unsigned char red   = pixel_cache[index].red;
            unsigned char green = pixel_cache[index].green;
            unsigned char blue  = pixel_cache[index].blue;

            ret( row, col ) = { red, green, blue } ;
    }

    return ret;
}


PixelImage pixelimage(1,1);

auto get_pixel_color = []( double x, double y, ColorComp cc ) -> double
{
    double ny = floor( y * pixelimage.getwidth()  );
    double nx = floor( x * pixelimage.getheight() );

    Pixel p00 = pixelimage( ny+0, nx+0 );
    Pixel p01 = pixelimage( ny+0, nx+1 );
    Pixel p10 = pixelimage( ny+1, nx+0 );
    Pixel p11 = pixelimage( ny+1, nx+1 );

    double lx = x - nx;
    double ly = y - ny;

    return lx * ly * p00(cc) + (1-lx) * ly * p01(cc) + lx * (1-ly) * p10(cc) + (1-lx) * (1-ly) * p11(cc);
};









int main() {
    
    // Initialize ImageMagick library
    InitializeMagick(nullptr);

    // Load the image file
    Image image;
    
    image.read("lena_color.tiff");

    // Get the dimensions of the image
    size_t width = image.columns();
    size_t height = image.rows();

    // Allocate stack-allocated arrays to store the pixel values
    unsigned char pixels_r[height][width];
    unsigned char pixels_g[height][width];
    unsigned char pixels_b[height][width];

    // Get the pixel values and store them in the arrays
    Pixels view(image);
    PixelPacket* pixel_cache = view.get(0, 0, width, height);
    for (size_t row = 0; row < height; row++) {
        for (size_t col = 0; col < width; col++) {
            size_t index = row * width + col;
            pixels_r[row][col] = pixel_cache[index].red;
            pixels_g[row][col] = pixel_cache[index].green;
            pixels_b[row][col] = pixel_cache[index].blue;

            pixel_cache[index].red   = 255 - pixel_cache[index].red;
            pixel_cache[index].green = 255 - pixel_cache[index].green;
            pixel_cache[index].blue  = 255 - pixel_cache[index].blue;
        }
    }

    view.sync();

    image.write("enal_color.tiff");

    // Print the pixel values
    for (size_t row = 0; row < height; row++) {
        for (size_t col = 0; col < width; col++) {
            std::cout << "(" << (int)pixels_r[row][col] << "," << (int)pixels_g[row][col] << "," << (int)pixels_b[row][col] << ") ";
        }
        std::cout << std::endl;
    }

    // Terminate ImageMagick library    
    return 0;
}
