
#include "vtkwriter.hpp"


VTK_MeshWriter::VTK_MeshWriter( SimplicialMesh& sm, std::ostream& os )
: mesh(sm), os(os)
{
    sm.check();
}


void VTK_MeshWriter::writePreamble( const char* name )
{
    // std::ostream& os = std::cout;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter::writeCoordinateBlock()
{
    // std::ostream& os = std::cout;
    os << "POINTS " << mesh.countsimplices(0) << " double" << nl;
    for( int v = 0; v < mesh.countsimplices(0); v++ )
        if( mesh.getouterdimension() == 2 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << mesh.getcoordinates().getdata(v,1) 
               << space 
               << 0.0
               << nl;
        } else if( mesh.getouterdimension() == 3 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << mesh.getcoordinates().getdata(v,1) 
               << space 
               << mesh.getcoordinates().getdata(v,2) 
               << nl;
        } else
            assert(false);
    
    os << nl;
}
        
        
void VTK_MeshWriter::writeTopDimensionalCells()
{
    const int innerdim = mesh.getinnerdimension();
    // std::ostream& os = std::cout;
    os << "CELLS " 
       << mesh.countsimplices(innerdim)
       << space 
       << innerdim+1 
       << nl;
    for( int S = 0; S < mesh.countsimplices(innerdim); S++ ) {
        for( int v = 0; v <= innerdim; v++ )
            os << mesh.getsubsimplices( innerdim, S, 0 )[v] << space;
        os << nl;
    }
    for( int S = 0; S < mesh.countsimplices(innerdim); S++ )
        os << ( innerdim == 3 ? 10 : 5 )
           << nl;
    
    os << nl;
}
        