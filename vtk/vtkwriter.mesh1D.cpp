
#include "vtkwriter.mesh1D.hpp"


VTK_MeshWriter_Mesh1D::VTK_MeshWriter_Mesh1D( Mesh& m1D, std::ostream& os )
: mesh(m1D), os(os)
{
    m1D.check();
    assert( m1D.getouterdimension() == 2 || m1D.getouterdimension() == 3 );
    assert( m1D.getinnerdimension() == 2 );
    assert( m1D.dimension_counted(1) );
    assert( m1D.dimension_counted(0) );
    assert( m1D.subsimplices_listed( 1, 0 ) );
        
}


void VTK_MeshWriter_Mesh1D::writePreamble( const char* name )
{
    // std::ostream& os = std::clog;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter_Mesh1D::writeCoordinateBlock()
{
    // std::ostream& os = std::clog;
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    for( int v = 0; v < mesh.count_simplices(0); v++ )
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
          unreachable();
        
    
    os << nl;
}
        
        
void VTK_MeshWriter_Mesh1D::writeTopDimensionalCells()
{
    // std::ostream& os = std::clog;
    
    os << "CELLS " 
       << mesh.count_simplices(1)
       << space 
       << 2 * mesh.count_simplices(1) + mesh.count_simplices(1)
       << nl;
    
    for( int S = 0; S < mesh.count_simplices(1); S++ ) 
        os << 2 
           << space
           << mesh.getsubsimplices(1,0,S)[0]
           << space
           << mesh.getsubsimplices(1,0,S)[1]
           << nl;
    
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.count_simplices(1) << nl;
    
    for( int S = 0; S < mesh.count_simplices(1); S++ )
        os << 3 << nl;
    
    os << nl;
}



void VTK_MeshWriter_Mesh1D::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.count_simplices(0) );
    
    os << "POINT_DATA " << mesh.count_simplices(0) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        os << scaling * data.at(v) << nl;
    
    os << std::endl;
    
}


void VTK_MeshWriter_Mesh1D::writeCellScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.count_simplices(1) );
    
    os << "CELL_DATA " << mesh.count_simplices(1) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // VECTORS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1) 
    os << "LOOKUP_TABLE default" << nl;
    
    for( int c = 0; c < mesh.count_simplices(1); c++ )
        os << scaling * data.at(c) << nl;
    
    os << std::endl;
    
}


void VTK_MeshWriter_Mesh1D::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_simplices(1) );
    
    os << "CELL_DATA " << mesh.count_simplices(1) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_simplices(1); c++ )
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    
    os << std::endl;
    
}





