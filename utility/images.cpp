#include <Magick++.h>
#include <iostream>

using namespace Magick;

// g++ images.cpp `Magick++-config --cxxflags --cppflags --ldflags --libs`
// sudo apt install graphicsmagick-libmagick-dev-compat libgraphicsmagick++1-dev libmagick++-6-headers
// sudo apt-get install libmagick++-dev
// sudo apt-get update; sudo apt-get upgrade; sudo apt-get install imagemagick-common

int main() {
    // Initialize ImageMagick library
    InitializeMagick(nullptr);

    // Load the BMP image file
    Image image("lena_color.tiff");

    // Get the dimensions of the image
    size_t width = image.columns();
    size_t height = image.rows();

    // Allocate stack-allocated arrays to store the pixel values
    unsigned char pixels[height][width][3];

    // Get the pixel values and store them in the arrays
    Pixels view(image);
    PixelPacket* pixel_cache = view.get(0, 0, width, height);
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            size_t index = y * width + x;
            pixels[y][x][0] = pixel_cache[index].red;
            pixels[y][x][1] = pixel_cache[index].green;
            pixels[y][x][2] = pixel_cache[index].blue;
        }
    }

    // Print the pixel values
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            std::cout << "(" << (int)pixels[y][x][0] << "," << (int)pixels[y][x][1] << "," << (int)pixels[y][x][2] << ") ";
        }
        std::cout << std::endl;
    }

    // Terminate ImageMagick library
    
    //TerminateMagick();

    return 0;
}
