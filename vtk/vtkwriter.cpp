
#include "vtkwriter.hpp"


VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name )
: VTKWriter( m, os, name, [&](int i) -> Float { return 0.; } )
{}

VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const FloatVector& z )
: VTKWriter( m, os, name, [&](int i) -> Float { return z[i]; } )
{}

VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const std::function<Float(int)>& func_z )
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

    writeCoordinateBlock( func_z );

    assert( current_stage == Stage::coordinate );

    writeTopDimensionalCells();
    
    assert( current_stage == Stage::cells );

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
    writeCoordinateBlock( [](int)->Float { return 0.; } );
    return *this;
}
        
VTKWriter VTKWriter::writeCoordinateBlock( const FloatVector& z )
{
    assert( mesh.count_simplices(0) == z.getdimension() );
    writeCoordinateBlock( [&](int i)->Float { return z.at(i); } );
    return *this;
}
        
        
VTKWriter VTKWriter::writeCoordinateBlock( const std::function<Float(int)>& func_z)
{
    assert( current_stage == Stage::preamble );
    current_stage = Stage::coordinate;
    
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        if( mesh.getouterdimension() == 1 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << 0.0 
               << space 
               << func_z(v)
               << nl;
        } else if( mesh.getouterdimension() == 2 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << mesh.getcoordinates().getdata(v,1) 
               << space 
               << func_z(v)
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
        
        
VTKWriter VTKWriter::writeTopDimensionalCells()
{
    assert( current_stage == Stage::coordinate);
    current_stage = Stage::cells;
    
    int topdim = mesh.getinnerdimension();
    
    assert( 1 <= topdim and topdim <= 3 );
    
    os << "CELLS " 
       << mesh.count_simplices(topdim)
       << space 
       << (topdim+1) * mesh.count_simplices(topdim) + mesh.count_simplices(topdim)
       << nl;
    
    for( int S = 0; S < mesh.count_simplices(topdim); S++ ) {
        os << topdim+1; 
        const auto list_of_vertices = mesh.getsubsimplices(topdim,0,S);
        for( int i = 0; i <= topdim; i++ ) os << space << list_of_vertices[i];
        os << nl;
    }
    
    os << nl;
    
    os << "CELL_TYPES" << space << mesh.count_simplices(topdim) << nl;
    
    if( topdim == 1 ) for( int S = 0; S < mesh.count_simplices(1); S++ ) os <<  3 << nl;
    if( topdim == 2 ) for( int S = 0; S < mesh.count_simplices(2); S++ ) os <<  5 << nl;
    if( topdim == 3 ) for( int S = 0; S < mesh.count_simplices(3); S++ ) os << 10 << nl;
    
    os << nl;

    return *this;
}

















VTKWriter VTKWriter::writeVertexScalarData( const std::function<Float(int)>& datafunction, const std::string name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::vertexdata );
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    if( current_stage != Stage::vertexdata ){
        os << "POINT_DATA " << mesh.count_simplices(0) << nl << nl;
        current_stage = Stage::vertexdata;
    }
    
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        os << scaling * datafunction(v) << nl;
    
    os << nl;

    return *this;
}

VTKWriter VTKWriter::writeCellScalarData( const std::function<Float(int)>& datafunction, const std::string name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::celldata );
    
    const int topdim = mesh.getinnerdimension();
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    if( current_stage != Stage::celldata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::celldata;
    }
    
    os << "SCALARS " << name << " double 1" << nl;
    // VECTORS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1) 
    os << "LOOKUP_TABLE default" << nl;
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ )
        os << scaling * datafunction(c) << nl;
    
    os << nl;

    return *this;
}
















VTKWriter VTKWriter::writeVertexScalarData( const std::function<Float(const FloatVector&)>& function, const std::string name, Float scaling )
{
    auto datafunction = [&](int c) -> Float { return function( mesh.get_midpoint(0,c) ); };

    writeVertexScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeVertexScalarData( const FloatVector& data, const std::string name, Float scaling )
{
    assert( data.getdimension() == mesh.count_simplices(0) );

    auto datafunction = [&](int c) -> Float { return data.at(c); };

    writeVertexScalarData( datafunction, name, scaling );
    
    return *this;
}






VTKWriter VTKWriter::writeCellScalarData( const std::function<Float(const FloatVector&)>& function, const std::string name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    auto datafunction = [&](int c) -> Float { return function( mesh.get_midpoint(topdim,c) ); };

    writeCellScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeCellScalarData( const FloatVector& data, const std::string name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    assert( data.getdimension() == mesh.count_simplices(topdim) );

    auto datafunction = [&](int c) -> Float { return data.at(c); };

    writeCellScalarData( datafunction, name, scaling );
    
    return *this;
}













VTKWriter VTKWriter::writeCellVectorData( const std::function<FloatVector(int)>& datafunction, const std::string name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::celldata );
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    const int topdim = mesh.getinnerdimension();
    
    if( current_stage != Stage::celldata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::celldata;
    }
    
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ ) {

        auto data_item = scaling * datafunction(c);

        os << data_item.at(0) << space 
           << ( mesh.getouterdimension() >= 2 ? data_item.at(1) : 0. ) << space 
           << ( mesh.getouterdimension() >= 3 ? data_item.at(2) : 0. ) << nl;
    }
    
    os << nl;

    return *this;
}


VTKWriter VTKWriter::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const std::string name, 
    Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_simplices(topdim) );
    
    auto datafunction = [&](int c) -> FloatVector { return FloatVector({datax[c],datay[c],dataz[c]}); };

    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}



VTKWriter VTKWriter::writeCellVectorData( const std::function<FloatVector(const FloatVector&)>& function, const std::string name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    auto datafunction = [&](int c) -> FloatVector { return function( mesh.get_midpoint(topdim,c) ); };
    
    writeCellVectorData( datafunction, name, scaling );

    return *this;
}















VTKWriter VTKWriter::writeCellVectorData_Whitney( const FloatVector& v, const std::string name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( v.getdimension() == mesh.count_simplices(topdim) * (mesh.getinnerdimension()+1) );
        
    auto datafunction = [&](int c) -> FloatVector {
    
        FloatVector coefficients( topdim+1 );
        
        for( int i = 0; i <= mesh.getinnerdimension(); i++ ) coefficients[i] = v.at( c * (topdim+1) + i );
        
        FloatVector directions = scaling * mesh.getGradientMatrix(topdim,c) * coefficients;

        return directions;
    };
    
    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}


VTKWriter VTKWriter::writeCellVectorData_Euclidean( int outerdim, const FloatVector& v, const std::string name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( v.getdimension() == mesh.count_simplices(topdim) * (mesh.getinnerdimension()+1) );
        
    auto datafunction = [&](int c) -> FloatVector {    
        return v.getslice( outerdim*c, outerdim );
    };
    
    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}




