
#include "vtkwriter.hpp"


VTKWriter::VTKWriter( Mesh& mesh, std::ostream& os, const std::string& name )
: mesh(mesh), os(os)
{
    mesh.check();
    
    assert( mesh.getouterdimension() == 1 || mesh.getouterdimension() == 2 || mesh.getouterdimension() == 3 );
    assert( mesh.getinnerdimension() == 1 || mesh.getinnerdimension() == 2 || mesh.getinnerdimension() == 3 );
    
    assert( mesh.dimension_counted( mesh.getinnerdimension() ) );
    assert( mesh.dimension_counted(0) );
    assert( mesh.subsimplices_listed( mesh.getinnerdimension(), 0 ) );
    
    writePreamble( name );
}


VTKWriter VTKWriter::writePreamble( const std::string& name )
{
    // std::ostream& os = std::clog;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
    
    return *this;
}



VTKWriter VTKWriter::writeCoordinateBlock()
{
    // std::ostream& os = std::clog;
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    for( int v = 0; v < mesh.count_simplices(0); v++ )
      if( mesh.getouterdimension() == 1 ) {
          os << mesh.getcoordinates().getdata(v,0)
              << space
              << 0.0 
              << space 
              << 0.0
              << nl;
      } else if( mesh.getouterdimension() == 2 ) {
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
      } else {
          unreachable();
      }
        
    
    os << nl;

    return *this;
}
        
        
VTKWriter VTKWriter::writeCoordinateBlock( const FloatVector& z )
{
    // std::ostream& os = std::clog;
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    assert( z.getdimension() == mesh.count_simplices(0) );
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
      if( mesh.getouterdimension() == 1 ) {
          os << mesh.getcoordinates().getdata(v,0)
              << space
              << 0.0 
              << space 
              << z.at(v)
              << nl;
      } else if( mesh.getouterdimension() == 2 ) {
          os << mesh.getcoordinates().getdata(v,0)
              << space
              << mesh.getcoordinates().getdata(v,1) 
              << space 
              << z.at(v)
              << nl;
      } else {
          unreachable();
      }
        
    
    os << nl;

    return *this;
}
        
        
VTKWriter VTKWriter::writeTopDimensionalCells()
{
    // std::ostream& os = std::clog;
    
    int topdim = mesh.getinnerdimension();
    
    assert( 1 <= topdim and topdim <= 3 );
    
    os << "CELLS " 
       << mesh.count_simplices(topdim)
       << space 
       << (topdim+1) * mesh.count_simplices(topdim) + mesh.count_simplices(topdim)
       << nl;
    
    for( int S = 0; S < mesh.count_simplices(topdim); S++ ) {
        os << topdim+1; 
        for( int i = 0; i <= topdim; i++ ) os << space << mesh.getsubsimplices(topdim,0,S)[i];
        os << nl;
    }
    
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.count_simplices(topdim) << nl;
    
    if( topdim == 1 ) for( int S = 0; S < mesh.count_simplices(1); S++ ) os <<  3 << nl;
    if( topdim == 2 ) for( int S = 0; S < mesh.count_simplices(2); S++ ) os <<  5 << nl;
    if( topdim == 3 ) for( int S = 0; S < mesh.count_simplices(3); S++ ) os << 10 << nl;
    
    os << nl;

    return *this;
}



VTKWriter VTKWriter::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
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

    return *this;
}


VTKWriter VTKWriter::writeCellScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    int topdim = mesh.getinnerdimension();
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.count_simplices(topdim) );
    
    os << "CELL_DATA " << mesh.count_simplices(topdim) << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // VECTORS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1) 
    os << "LOOKUP_TABLE default" << nl;
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ )
        os << scaling * data.at(c) << nl;
    
    os << std::endl;

    return *this;
}


VTKWriter VTKWriter::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    int topdim = mesh.getinnerdimension();
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_simplices(topdim) );
    
    os << "CELL_DATA " << mesh.count_simplices(topdim) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ )
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    
    os << std::endl;

    return *this;
}





