
#include "vtkwriter.manifold2D.hpp"


VTK_MeshWriter_Manifold2D::VTK_MeshWriter_Manifold2D( MeshManifold2D& m2d, std::ostream& os )
: mesh(m2d), os(os)
{
    m2d.check();
}


void VTK_MeshWriter_Manifold2D::writePreamble( const char* name )
{
    // std::ostream& os = std::cout;
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
}


void VTK_MeshWriter_Manifold2D::writeCoordinateBlock()
{
    // std::ostream& os = std::cout;
    os << "POINTS " << mesh.count_vertices() << " double" << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ )
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
        
        
void VTK_MeshWriter_Manifold2D::writeTopDimensionalCells()
{
    const int innerdim = 2;
    // std::ostream& os = std::cout;
    
    os << "CELLS " 
       << mesh.count_triangles()
       << space 
       << 3 * mesh.count_triangles() + mesh.count_triangles()
       << nl;

    for( int S = 0; S < mesh.count_triangles(); S++ ) 
        os << 3 
           << space
           << mesh.get_triangle_vertices(S)[0]
           << space
           << mesh.get_triangle_vertices(S)[1]
           << space
           << mesh.get_triangle_vertices(S)[2]
           << nl;
    
    os << std::endl;
    
    os << "CELL_TYPES" << space << mesh.count_triangles() << nl;
    
    for( int S = 0; S < mesh.count_triangles(); S++ )
        os << ( innerdim == 3 ? 10 : 5 )
           << nl;
    
    os << nl;
}



void VTK_MeshWriter_Manifold2D::writeVertexScalarData( const FloatVector& data, const char* name, Float scaling )
{
    
    assert( name != nullptr );
    assert( data.getdimension() == mesh.count_vertices() );
    
    os << "POINT_DATA " << mesh.count_vertices() << nl;
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.count_vertices(); v++ )
        os << scaling * data.at(v) << nl;
    
    os << std::endl;
    
}


void VTK_MeshWriter_Manifold2D::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const char* name, 
    Float scaling )
{
    
    assert( name != nullptr );
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_triangles() );
    
    const int innerdim = 2;
    
    os << "CELL_DATA " << mesh.count_triangles() << nl;
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_triangles(); c++ )
    {
      os << scaling * datax.at(c) << space 
         << scaling * datay.at(c) << space 
         << scaling * dataz.at(c) 
         << nl;
    }
    
    os << std::endl;
    
}





