
#include "vtkwriter.mesh3D.hpp"


VTK_MeshWriter_Mesh3D::VTK_MeshWriter_Mesh3D( Mesh& m3D, std::ostream& os )
: mesh(m3D), os(os)
{
    m3D.check();
    assert( m3D.getouterdimension() == 3 );
    assert( m3D.getinnerdimension() == 3 );
    assert( m3D.dimensioncounted(3) );
    assert( m3D.dimensioncounted(0) );
    assert( m3D.subsimplices_listed( 3, 0 ) );
}


void VTK_MeshWriter_Mesh3D::writePreamble( const char* name )
{
    // std::ostream& os = std::cout;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter_Mesh3D::writeCoordinateBlock()
{
    // std::ostream& os = std::cout;
    os << "POINTS " << mesh.countsimplices(0) << " double" << nl;
    for( int v = 0; v < mesh.countsimplices(0); v++ )
      os << mesh.getcoordinates().getdata(v,0)
         << space
         << mesh.getcoordinates().getdata(v,1) 
         << space 
         << mesh.getcoordinates().getdata(v,2) 
         << nl;
      
    os << nl;
}
        
        
void VTK_MeshWriter_Mesh3D::writeTopDimensionalCells()
{
    // std::ostream& os = std::cout;
    
    os << "CELLS " 
       << mesh.countsimplices(3)
       << space 
       << 4 * mesh.countsimplices(3) + mesh.countsimplices(3)
       << nl;
    
    for( int S = 0; S < mesh.countsimplices(3); S++ ) 
        os << 4 
           << space
           << mesh.getsubsimplices(3,0,S)[0]
           << space
           << mesh.getsubsimplices(3,0,S)[1]
           << space
           << mesh.getsubsimplices(3,0,S)[2]
           << space
           << mesh.getsubsimplices(3,0,S)[3]
           << nl;
    
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.countsimplices(3) << nl;
    
    for( int S = 0; S < mesh.countsimplices(3); S++ )
        os << 10 << nl;
    
    os << nl;
}



void VTK_MeshWriter_Mesh3D::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.countsimplices(0) );
    
    os << "POINT_DATA " << mesh.countsimplices(0) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.countsimplices(0); v++ )
        os << scaling * data.at(v) << nl;
    
    os << std::endl;
    
}


void VTK_MeshWriter_Mesh3D::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.countsimplices(3) );
    
    os << "CELL_DATA " << mesh.countsimplices(3) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.countsimplices(3); c++ )
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    
    os << std::endl;
    
}





