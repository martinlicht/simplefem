

/**/

#include <iostream>
#include <fstream>

#include "../../matrixmarket/matrixmarket.hpp"


using namespace std;

int main()
{
    cout << "Unit Test for Matrix Market Reading" << endl;
    
    {
        
        int rows, columns;
        std::vector<MatrixMarket::MatrixMarketEntry> entries;
        
        ifstream fs( "gr_30_30.mtx" );
        
        assert( fs.good() );
        
        MatrixMarket::Read( fs, rows, columns, entries );
        
        std::cout << rows << ' ' << columns << ' ' << std::endl;
        for( MatrixMarket::MatrixMarketEntry mme : entries )
            std::cout << mme.row << ' ' << mme.column << ' ' << mme.value << std::endl;
        
        fs.close();
        
//         // MeshSimplicial1D M = UnitCubeTriangulation(3,3);
//         MeshSimplicial1D M = StandardInterval1D();
//         
//         cout << M << endl;
//         
//         cout << "Print VTK-type file" << endl;
//         
//         
//         {
//             
//             VTK_MeshWriter_Mesh1D vtk( M, fs );
//             vtk.writePreamble( "Mein erster Test" );
//             vtk.writeCoordinateBlock();
//             vtk.writeTopDimensionalCells();
//             
//         }
//         
//      
    }
    
    cout << "Finished Unit Test" << endl;
    
    return 0;
}
