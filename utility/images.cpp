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

std::string render_number( double num, int tail )
{
    int lead = leading_digits( num );
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
    <polygon points="10,20 30,40 50,60" fill="#FF3399" stroke="none"></polygon>
    <polygon points="10,20 40,50 0010,10" fill="#FF3399" stroke="#11CCDD" stroke-width="1"></polygon>

    #include <iostream>
    #include <sstream>
    #include <string>

    std::string generatePolygonString(
        double x1, double y1, double x2, double y2, double x3, double y3, 
        std::string fill, 
        std::string stroke,
        double stroke_width
    ) {
        std::stringstream ss;
        ss << "<polygon points=\"" 
           << x1 << "," << y1 << " " 
           << x2 << "," << y2 << " " 
           << x3 << "," << y3 << "\""  
        ss << " fill=\"" << fill << "\""  
        ss << " stroke=\"" << stroke << "\""  
        ss << " stroke-width=\"" << stroke_width << "\""  
        ss << "></polygon>";
        return ss.str();
    }

    int main() {
        std::string polygonString = generatePolygonString(10, 20, 30, 40, 50, 60, "#FF3399", "none");
        std::cout << polygonString << std::endl;
        return 0;
    }


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

enum class ColorChannel { red, green, blue };

struct PixelColor
{
    unsigned char red;
    unsigned char green;
    unsigned char blue;

    unsigned char operator()(ColorChannel cc) const { 
        switch(cc) {
            case ColorChannel::red:   return red;
            case ColorChannel::green: return green;
            case ColorChannel::blue:  return blue;
        }
    }
};

class PixelImage
{
    private:
        size_t height;
        size_t width;
        std::vector<PixelColor> data;

    public:

        PixelImage( size_t height, size_t width ) 
        : height(height), width(width), data( width*height, {128,128,128} ) 
        {}
        
        ~PixelImage() {}

        size_t getheight() const { return height; };

        size_t getwidth() const { return width; };

        PixelColor& operator()( size_t row, size_t col ) {
            assert( row >= 0 && row < height && col >= 0 && col < width );
            return data.at( row * width + col );
        }
        
        PixelColor operator()( size_t row, size_t col ) const {
            assert( row >= 0 && row < height && col >= 0 && col < width );
            return data.at( row * width + col );
        }

        auto get_interpolated_function() const {
            return [&]( double x, double y, ColorChannel cc ) -> double
            {
                const PixelImage& pixelimage = *this;
                
                double ny = floor( y * pixelimage.getwidth()  );
                double nx = floor( x * pixelimage.getheight() );
                double lx = x - nx;
                double ly = y - ny;

                PixelColor p00 = pixelimage( ny+0, nx+0 );
                PixelColor p01 = pixelimage( ny+0, nx+1 );
                PixelColor p10 = pixelimage( ny+1, nx+0 );
                PixelColor p11 = pixelimage( ny+1, nx+1 );

                return lx * ly * p00(cc) + (1-lx) * ly * p01(cc) + lx * (1-ly) * p10(cc) + (1-lx) * (1-ly) * p11(cc);
            };
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
    
    for( size_t row = 0; row < height; row++) 
    for( size_t col = 0; col < width; col++)
    {
        ColorRGB color = image.pixelColor( row, col ); 
        ret(row,col).red   = (int)(255*color.red());
        ret(row,col).green = (int)(255*color.green());
        ret(row,col).blue  = (int)(255*color.blue());
    }

    return ret;

    // Pixels view(image);
    // const PixelPacket* pixel_cache = view.getConst(0, 0, width, height);
    // for( size_t row = 0; row < height; row++) 
    // for( size_t col = 0; col < width; col++) 
    // {
    //         size_t index = row * width + col;

    //         unsigned char red   = pixel_cache[index].red;
    //         unsigned char green = pixel_cache[index].green;
    //         unsigned char blue  = pixel_cache[index].blue;

    //         ret( row, col ) = { red, green, blue } ;
    // }

    // return ret;
}

void savePixelImage( const PixelImage& pim, std::string str )
{
    // InitializeMagick(nullptr);

    size_t width  = pim.getwidth();
    size_t height = pim.getheight();

    Image image( Geometry(width, height), "white" );

    for( size_t row = 0; row < height; row++) 
    for( size_t col = 0; col < width; col++)
    {
        auto color = pim( row, col );
        image.pixelColor( row, col, ColorRGB( color.red/255., color.green/255., color.blue/255. ) ); 
    }
        
    image.write( str ); 
    return;
}








int main() {
    
    auto pim = readPixelImage("lena_color.tiff");

    for( size_t row = 0; row < pim.getheight(); row++) 
    for( size_t col = 0; col < pim.getwidth(); col++) 
    {
        pim(row,col).red = 255 - pim(row,col).red;
        pim(row,col).green = 255 - pim(row,col).green;
        pim(row,col).blue = 255 - pim(row,col).blue;
    }

    savePixelImage(pim,"anel_color.tiff");

    return 0;

    // // Initialize ImageMagick library
    // InitializeMagick(nullptr);

    // // Load the image file
    // Image image;
    
    // image.read("lena_color.tiff");

    // // Get the dimensions of the image
    // size_t width = image.columns();
    // size_t height = image.rows();

    // // Allocate stack-allocated arrays to store the pixel values
    // unsigned char pixels_r[height][width];
    // unsigned char pixels_g[height][width];
    // unsigned char pixels_b[height][width];

    // // Get the pixel values and store them in the arrays
    // Pixels view(image);
    // PixelPacket* pixel_cache = view.get(0, 0, width, height);
    // for( size_t row = 0; row < height; row++) {
    //     for( size_t col = 0; col < width; col++) {
    //         size_t index = row * width + col;
    //         pixels_r[row][col] = pixel_cache[index].red;
    //         pixels_g[row][col] = pixel_cache[index].green;
    //         pixels_b[row][col] = pixel_cache[index].blue;

    //         pixel_cache[index].red   = 255 - pixel_cache[index].red;
    //         pixel_cache[index].green = 255 - pixel_cache[index].green;
    //         pixel_cache[index].blue  = 255 - pixel_cache[index].blue;
    //     }
    // }

    // view.sync();

    // image.write("enal_color.tiff");

    // // Print the pixel values
    // for( size_t row = 0; row < height; row++) {
    //     for( size_t col = 0; col < width; col++) {
    //         std::cout << "(" << (int)pixels_r[row][col] << "," << (int)pixels_g[row][col] << "," << (int)pixels_b[row][col] << ") ";
    //     }
    //     std::cout << std::endl;
    // }

    // // Terminate ImageMagick library    
    // return 0;
}
