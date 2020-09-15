
#include "vtkwriter.hpp"


VTKWriter::VTKWriter( const Mesh& mesh, std::ostream& os, const std::string& name )
: mesh(mesh), os(os), current_stage(Stage::nothing)
{
    mesh.check();
    
    assert( mesh.getouterdimension() == 1 || mesh.getouterdimension() == 2 || mesh.getouterdimension() == 3 );
    assert( mesh.getinnerdimension() == 1 || mesh.getinnerdimension() == 2 || mesh.getinnerdimension() == 3 );
    
    assert( mesh.dimension_counted( mesh.getinnerdimension() ) );
    assert( mesh.dimension_counted(0) );
    assert( mesh.subsimplices_listed( mesh.getinnerdimension(), 0 ) );
    
    writePreamble( name );
    
    assert( current_stage == Stage::preamble );
}


VTKWriter VTKWriter::writePreamble( const std::string& name )
{
    assert( current_stage == Stage::nothing );
    current_stage = Stage::preamble;
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    // assert( count_white_space( name ) == 0 );
    
    
    // std::ostream& os = std::clog;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
    
    return *this;
}



VTKWriter VTKWriter::writeCoordinateBlock()
{
    assert( current_stage == Stage::preamble );
    current_stage = Stage::coordinate;
    
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
    assert( current_stage == Stage::preamble );
    current_stage = Stage::coordinate;
    
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
    assert( current_stage == Stage::coordinate);
    current_stage = Stage::cells;
    
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



VTKWriter VTKWriter::writeVertexScalarData( const FloatVector& data, const std::string name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::vertexdata );
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    assert( data.getdimension() == mesh.count_simplices(0) );
    
    if( current_stage != Stage::vertexdata ){
        os << "POINT_DATA " << mesh.count_simplices(0) << nl << nl;
        current_stage = Stage::vertexdata;
    }
    
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        os << scaling * data.at(v) << nl;
    
    os << std::endl;

    return *this;
}


VTKWriter VTKWriter::writeCellScalarData( const FloatVector& data, const std::string name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::celldata );
    
    const int topdim = mesh.getinnerdimension();
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    assert( data.getdimension() == mesh.count_simplices(topdim) );
    
    if( current_stage != Stage::celldata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::celldata;
    }
    
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
    const std::string name, 
    Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::celldata );
    
    const int topdim = mesh.getinnerdimension();
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_simplices(topdim) );
    
    if( current_stage != Stage::celldata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::celldata;
    }
    
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





