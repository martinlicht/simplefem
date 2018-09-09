
#include "vtkwriter.hpp"


VTK_MeshWriter::VTK_MeshWriter( SimplicialMesh& sm, std::ostream& os )
: mesh(sm), os(os)
{
    sm.check();
}


void VTK_MeshWriter::writePreamble( const char* name )
{
    // std::ostream& os = std::clog;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter::writeCoordinateBlock()
{
    // std::ostream& os = std::clog;
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
    // std::ostream& os = std::clog;
    
    os << "CELLS " 
       << mesh.countsimplices(innerdim)
       << space 
       << ( innerdim+1 ) * mesh.countsimplices(innerdim) + mesh.countsimplices(innerdim)
       << nl;
    for( int S = 0; S < mesh.countsimplices(innerdim); S++ ) {
        os << innerdim+1 
           << space;
        for( int v = 0; v <= innerdim; v++ )
            os << mesh.getsubsimplices( innerdim, S, 0 )[v]
               << space;
        os << nl;
    }
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.countsimplices(innerdim) << nl;
    for( int S = 0; S < mesh.countsimplices(innerdim); S++ )
        os << ( innerdim == 3 ? 10 : 5 )
           << nl;
    
    os << nl;
}



void VTK_MeshWriter::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.countsimplices(0) );
    
    os << "POINT_DATA " << mesh.countsimplices(0) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.countsimplices(0); v++ )
    {
      os << scaling * data.at(v) << nl;
    }
    
    os << std::endl;
    
}


void VTK_MeshWriter::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.countsimplices( mesh.getinnerdimension() ) );
    
    const int innerdim = mesh.getinnerdimension();
    
    os << "CELL_DATA " << mesh.countsimplices( innerdim ) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.countsimplices( innerdim ); c++ )
    {
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    }
    
    os << std::endl;
    
}





