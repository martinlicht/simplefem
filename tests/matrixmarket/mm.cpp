

/**/

#include <iostream>
#include <fstream>

#include "../../matrixmarket/matrixmarket.hpp"


using namespace std;

int main()
{
    LOG << "Unit Test for Matrix Market Reading" << endl;
    
    {
        
        int rows, columns;
        std::vector<MatrixMarket::MatrixMarketEntry> entries;
        
        ifstream fs( "gr_30_30.mtx" );
        
        assert( fs.good() );
        
        MatrixMarket::Read( fs, rows, columns, entries );
        
        LOG << rows << ' ' << columns << ' ' << std::endl;
        for( MatrixMarket::MatrixMarketEntry mme : entries )
            LOG << mme.row << ' ' << mme.column << ' ' << mme.value << std::endl;
        
        fs.close();
        
//         // MeshSimplicial1D M = UnitCubeTriangulation(3,3);
//         MeshSimplicial1D M = StandardInterval1D();
//         
//         LOG << M << endl;
//         
//         LOG << "Print VTK-type file" << endl;
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
    
    LOG << "Finished Unit Test" << endl;
    
    return 0;
}
