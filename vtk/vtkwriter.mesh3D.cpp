
#include "vtkwriter.mesh3D.hpp"


VTK_MeshWriter_Mesh3D::VTK_MeshWriter_Mesh3D( Mesh& m3D, std::ostream& os )
: mesh(m3D), os(os)
{
    m3D.check();
    assert( m3D.getouterdimension() == 3 );
    assert( m3D.getinnerdimension() == 3 );
    assert( m3D.dimension_counted(3) );
    assert( m3D.dimension_counted(0) );
    assert( m3D.subsimplices_listed( 3, 0 ) );
}


void VTK_MeshWriter_Mesh3D::writePreamble( const std::string& str )
{
    writePreamble( str.c_str() );
}

void VTK_MeshWriter_Mesh3D::writePreamble( const char* name )
{
    // std::ostream& os = std::clog;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter_Mesh3D::writeCoordinateBlock()
{
    // std::ostream& os = std::clog;
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    for( int v = 0; v < mesh.count_simplices(0); v++ )
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
    // std::ostream& os = std::clog;
    
    os << "CELLS " 
       << mesh.count_simplices(3)
       << space 
       << 4 * mesh.count_simplices(3) + mesh.count_simplices(3)
       << nl;
    
    for( int S = 0; S < mesh.count_simplices(3); S++ ) 
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
    
    os << "CELL_TYPES" << space << mesh.count_simplices(3) << nl;
    
    for( int S = 0; S < mesh.count_simplices(3); S++ )
        os << 10 << nl;
    
    os << nl;
}



void VTK_MeshWriter_Mesh3D::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
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


void VTK_MeshWriter_Mesh3D::writeCellScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.count_simplices(3) );
    
    os << "CELL_DATA " << mesh.count_simplices(3) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // VECTORS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1) 
    os << "LOOKUP_TABLE default" << nl;
    
    for( int c = 0; c < mesh.count_simplices(3); c++ )
        os << scaling * data.at(c) << nl;
    
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
    assert( datax.getdimension() == mesh.count_simplices(3) );
    
    os << "CELL_DATA " << mesh.count_simplices(3) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_simplices(3); c++ )
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    
    os << std::endl;
    
}





