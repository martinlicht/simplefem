#ifndef INCLUDEGUARD_DENSE_MATRIX
#define INCLUDEGUARD_DENSE_MATRIX

#include <vector>

class DenseMatrix;

#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/scalingoperator.hpp"
#include "../operators/diagonaloperator.hpp"
#include "../sparse/sparsematrix.hpp"
// #include "matrixalgorithm.hpp"



/************************
****
****  Class for Dense Matrices 
****  - instantiates LinearOperator
****  - only basic linear arithmetics, and also matrix multiplication 
****  
************************/

class DenseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

    public:
        
        DenseMatrix( int dim, Float value = 0. );
        DenseMatrix( int dim, const std::function<Float(int,int)>& generator );
        DenseMatrix( int dim, const std::vector<FloatVector>& coldata );
        
        DenseMatrix( int rows, int columns, Float value = 0. );
        DenseMatrix( int rows, int columns, const std::function<Float(int,int)>& generator );
        DenseMatrix( int rows, int columns, const std::vector<FloatVector>& coldata );
        
        explicit DenseMatrix( const ScalingOperator& );
        explicit DenseMatrix( const DiagonalOperator& );
        explicit DenseMatrix( const SparseMatrix& );
        explicit DenseMatrix( const FloatVector& );
        
        virtual ~DenseMatrix();
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;
        
        virtual FloatVector apply( const FloatVector& add, Float scaling ) const override;
        
        /* matrix level point of view */
        
        int numrows() const;
        int numcolumns() const;
        
        /* Access entries */
        
        Float get(int,int) const;
        void set(int,int,Float);
        Float& at( int, int );
        const Float& at( int, int ) const;
        Float& operator()( int, int );
        const Float& operator()( int, int ) const;
        
        /* Access rows and columns */
        
        FloatVector getrow( int r ) const;
        FloatVector getcolumn( int c ) const;
        void setrow( int r, const FloatVector& row );
        void setcolumn( int c, const FloatVector& column );
        void addrow( int r, const FloatVector& row, Float s );
        void addcolumn( int c, const FloatVector& column, Float s );
        
        void swaprow( int r1, int r2 );
        void swapcolumn( int c1, int c2 );
        void scalerow( int r, Float alpha );
        void scalecolumn( int c, Float alpha );
        void addrow( int r1, int r2, Float alpha );
        void addcolumn( int c1, int c2, Float alpha );
        
        /* Flatten (and raise?) */
        
        FloatVector flattencolumns() const;
        FloatVector flattenrows() const;
        
        
        /* Generate standard matrices */
        
        void zeromatrix();
        void randommatrix();
        void randomintegermatrix( int min, int max );
        void unitmatrix();
        void indexmapping( const IndexMap& );
        
        /* Basic manipulation */
        
        void scale( Float );
        void set( Float );
        void add( Float );
        
        /* Special operations */
        
        DenseMatrix submatrix( const IndexMap& rows, const IndexMap& columns ) const;
        
        /* Arithmetic operations */
        
        void add( const DenseMatrix& );
        void add( Float, const DenseMatrix& );
        void add( Float, Float, const DenseMatrix& );
        
        /* Measurements */
        
        Float maxabsoluteentry() const;
        
        Float norm() const;
        
        Float maxnorm() const;
        
        Float sumnorm() const;
        
        Float lpnorm( Float ) const;
        
        
        /* Investigations */
        
        bool issquare() const;
        
        bool issymmetric() const;
        
        bool isantisymmetric() const;
        
        bool isfinite() const;
        
        bool iszero() const;
        
        bool ispositive() const;
        
        bool isnegative() const;
        
        bool isnonnegative() const;
        
        bool isnonpositive() const;
        
        
        
    private:
        
        std::vector<Float> entries;
        
};
  

inline DenseMatrix MatrixMult( const DenseMatrix& left, const DenseMatrix& right )
{
  left.check();
  right.check();
  
  const int lin = left.getdimin();
  const int lout = left.getdimout();
  const int rin = right.getdimin();
  const int rout = right.getdimout();
  
  assert( lin == rout );
  
  DenseMatrix ret( lout, rin );
  ret.zeromatrix();
  
  for( int lo = 0; lo < lout; lo++ )
  for( int ri = 0; ri < rin; ri++ )
  for( int m = 0; m < rout; m++ )
  ret( lo, ri ) += left( lo, m ) * right( m, ri );
  
  ret.check();
  return ret;
  
}


Float Determinant( const DenseMatrix& );

inline Float determinant( const DenseMatrix& mat ) 
{
  return Determinant(mat);
}





inline DenseMatrix& operator+=( DenseMatrix& left, const DenseMatrix& right )
{
    left.add( right );
    return left;
}

inline DenseMatrix& operator-=( DenseMatrix& left, const DenseMatrix& right )
{
    left.add( -1., right );
    return left;
}

inline DenseMatrix& operator*=( DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix temp( left );
    temp = MatrixMult( left, right );
    left = temp;
    return left;
}

inline DenseMatrix& operator*=( DenseMatrix& left, Float right )
{
    left.scale( right );
    return left;
}

inline DenseMatrix& operator/=( DenseMatrix& left, Float right )
{
    left.scale( 1. / right );
    return left;
}




inline DenseMatrix operator+( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret += right;
    return ret;
}

inline DenseMatrix operator-( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret -= right;
    return ret;
}

inline DenseMatrix operator*( const DenseMatrix& left, const DenseMatrix& right )
{
    DenseMatrix ret( left );
    ret *= right;
    return ret;
}

inline DenseMatrix operator*( Float left, const DenseMatrix& right )
{
    DenseMatrix ret( right );
    ret *= left;
    return ret;
}

inline DenseMatrix operator*( const DenseMatrix& left, Float right )
{
    DenseMatrix ret( left );
    ret *= right;
    return ret;
}

inline DenseMatrix operator/( const DenseMatrix& left, Float right )
{
    DenseMatrix ret( left );
    ret /= right;
    return ret;
}



#endif
