
#include "vtkwriter.manifold2D.hpp"


VTK_MeshWriter_Mesh2D::VTK_MeshWriter_Mesh2D( Mesh& m2d, std::ostream& os )
: mesh(m2d), os(os)
{
    m2d.check();
    assert( m2d.getouterdimension() == 2 || m2d.getouterdimension() == 3 );
    assert( m2d.getinnerdimension() == 2 );
    assert( m2d.dimensioncounted(2) );
    assert( m2d.dimensioncounted(0) );
    assert( m2d.subsimplices_listed( 2, 0 ) );
        
}


void VTK_MeshWriter_Mesh2D::writePreamble( const char* name )
{
    // std::ostream& os = std::cout;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter_Mesh2D::writeCoordinateBlock()
{
    // std::ostream& os = std::cout;
    os << "POINTS " << mesh.countsimplices(0) << " double" << nl;
    for( int v = 0; v < mesh.countsimplices(0); v++ )
    {
        
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
        
      
    }
    os << nl;
}
        
        
void VTK_MeshWriter_Mesh2D::writeTopDimensionalCells()
{
    const int innerdim = 2;
    // std::ostream& os = std::cout;
    
    os << "CELLS " 
       << mesh.countsimplices(2)
       << space 
       << 3 * mesh.countsimplices(2) + mesh.countsimplices(2)
       << nl;

    for( int S = 0; S < mesh.countsimplices(2); S++ ) 
        os << 3 
           << space
           << mesh.getsubsimplices(2,0,S)[0]
           << space
           << mesh.getsubsimplices(2,0,S)[1]
           << space
           << mesh.getsubsimplices(2,0,S)[2]
           << nl;
    
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.countsimplices(2) << nl;
    
    for( int S = 0; S < mesh.countsimplices(2); S++ )
        os << ( innerdim == 3 ? 10 : 5 )
           << nl;
    
    os << nl;
}



void VTK_MeshWriter_Mesh2D::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
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


void VTK_MeshWriter_Mesh2D::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.countsimplices(2) );
    
    const int innerdim = 2;
    
    os << "CELL_DATA " << mesh.countsimplices(2) << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.countsimplices(2); c++ )
    {
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    }
    
    os << std::endl;
    
}





