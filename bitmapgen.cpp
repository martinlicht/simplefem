#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cmath>

#include <array>
#include <vector>

#pragma pack(push, 1) // Ensure no padding in the structures
struct BMPFileHeader {
    uint16_t fileType = 0x4D42; // "BM"
    uint32_t fileSize;
    uint16_t reserved1 = 0;
    uint16_t reserved2 = 0;
    uint32_t offsetData;
};

struct BMPInfoHeader {
    uint32_t size = 40; // Header size in bytes
    int32_t width;
    int32_t height;
    uint16_t planes = 1;
    uint16_t bitCount = 24; // Bits per pixel
    uint32_t compression = 0;
    uint32_t sizeImage = 0;
    int32_t xPixelsPerMeter = 0;
    int32_t yPixelsPerMeter = 0;
    uint32_t colorsUsed = 0;
    uint32_t colorsImportant = 0;
};
#pragma pack(pop)

void writeBMP(const std::string& filename, int width, int height, std::vector<uint8_t>& pixelData) {
    // Define the bitmap headers
    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;

    infoHeader.width = width;
    infoHeader.height = height;

    int rowStride = (width * 3 + 3) & ~3; // Align rows to 4 bytes
    infoHeader.sizeImage = rowStride * height;
    fileHeader.fileSize = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader) + infoHeader.sizeImage;
    fileHeader.offsetData = sizeof(BMPFileHeader) + sizeof(BMPInfoHeader);

    // Open the file for writing
    std::ofstream outFile(filename, std::ios::binary);
    if (!outFile) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // Write headers
    outFile.write(reinterpret_cast<const char*>(&fileHeader), sizeof(fileHeader));
    outFile.write(reinterpret_cast<const char*>(&infoHeader), sizeof(infoHeader));

    // Write pixel data
    for (int y = 0; y < height; ++y) {
        outFile.write(reinterpret_cast<const char*>(&pixelData[y * rowStride]), rowStride);
    }

    outFile.close();
    std::cout << "BMP file written to " << filename << std::endl;
}

void drawLine(std::vector<uint8_t>& pixelData, int width, int height, int x1, int y1, int x2, int y2) {
    int rowStride = (width * 3 + 3) & ~3;

    // Bresenham's line algorithm
    int dx = std::abs(x2 - x1);
    int dy = std::abs(y2 - y1);
    int sx = (x1 < x2) ? 1 : -1;
    int sy = (y1 < y2) ? 1 : -1;
    int err = dx - dy;

    while (true) {
        if (x1 >= 0 && x1 < width && y1 >= 0 && y1 < height) {
            int offset = (y1 * rowStride) + (x1 * 3);
            pixelData[offset] = 255;     // Blue
            pixelData[offset + 1] = 255; // Green
            pixelData[offset + 2] = 255; // Red
        }

        if (x1 == x2 && y1 == y2) break;

        int e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x1 += sx;
        }
        if (e2 < dx) {
            err += dx;
            y1 += sy;
        }
    }
}

int main() {
    
    std::vector<double> coordinates 
    =
    {
        -1.0, -1.0, // 0
        -1.0,  1.0, // 1
            1.0, -1.0, // 2
            1.0,  1.0, // 3
            //
        -1.0, 0.1, // 4
        -0.2,-1.0, // 5
            0.3, 1.0,  // 6
            1.0,-0.2,  // 7
            //
        -0.3, -0.5, // 8
        -0.4,  0.4, // 9
            0.4, -0.5, // A
            0.5,  0.3, // B
            //
        
        };
    for( auto& coord : coordinates ) coord = ( coord + 1. ) / 2.;


    // Define triangles
    std::vector<std::array<int, 3>> triangles = {
        {0, 4, 8},
        {0, 5, 8},
        {5, 8, 10},
        {2, 5, 10},
        {2, 7, 10},
        {4, 8, 9},
        {8, 9, 11},
        {8, 10, 11},
        {7, 10, 11},
        {1, 4, 9},
        {1, 6, 9},
        {6, 9, 11},
        {3, 6, 11},
        {3, 7, 11}
    };
    
    int width = 1200;  // Hardcoded width
    int height = 1200; // Hardcoded height

    int rowStride = (width * 3 + 3) & ~3;
    std::vector<uint8_t> pixelData(rowStride * height, 0); // Blank canvas (black)

    for( const auto& t : triangles )
    {
        int v0 = t[0];
        int v1 = t[1];
        int v2 = t[2];

        double x0 = coordinates[ 2 * v0 + 0 ];
        double y0 = coordinates[ 2 * v0 + 1 ];
        
        double x1 = coordinates[ 2 * v1 + 0 ];
        double y1 = coordinates[ 2 * v1 + 1 ];
        
        double x2 = coordinates[ 2 * v2 + 0 ];
        double y2 = coordinates[ 2 * v2 + 1 ];

        int cx0 = 100 + 1000 * x0;
        int cy0 = 100 + 1000 * y0;
        
        int cx1 = 100 + 1000 * x1;
        int cy1 = 100 + 1000 * y1;
        
        int cx2 = 100 + 1000 * x2;
        int cy2 = 100 + 1000 * y2;
        
        drawLine(pixelData, width, height, cx0, cy0, cx1, cy1 );
        drawLine(pixelData, width, height, cx0, cy0, cx2, cy2 );
        drawLine(pixelData, width, height, cx1, cy1, cx2, cy2 );
        
    }

    
    
    
    
    
    // Draw some lines
    // drawLine(pixelData, width, height, 100, 100, 700, 500);
    // drawLine(pixelData, width, height, 400, 100, 400, 500);

    writeBMP("temp.bmp", width, height, pixelData);

    return 0;
}
