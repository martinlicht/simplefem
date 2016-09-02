#ifndef INCLUDEGUARD_DENSEMATRIX
#define INCLUDEGUARD_DENSEMATRIX

#include <vector>

class DenseMatrix;

#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "linearoperator.hpp"
#include "scalingoperator.hpp"
#include "diagonaloperator.hpp"
#include "sparsematrix.hpp"



/************************
****
****  Class for Dense Matrices 
****  - instantiates LinearOperator
****  - only basic linear arithmetics 
****  
************************/

class DenseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

    public:
        
        explicit DenseMatrix( int dim );
        DenseMatrix( int dim, std::function<Float(int,int)> generator );
        DenseMatrix( int rows, int columns );
        DenseMatrix( int rows, int columns, std::function<Float(int,int)> generator );
        explicit DenseMatrix( const ScalingOperator& );
        explicit DenseMatrix( const DiagonalOperator& );
        explicit DenseMatrix( const SparseMatrix& );
        explicit DenseMatrix( const FloatVector& );
        virtual ~DenseMatrix();
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        
        virtual void applyadd( FloatVector& dest, const FloatVector& add, Float s, Float t ) const override;
        
        /* Questions */
        
        bool issquare() const;
        
        /* Access entries */
        
        Float get(int,int) const;
        void set(int,int,Float);
        Float& operator()( int, int );
        const Float& operator()( int, int ) const;
        
        /* Access rows and columns */
        
        FloatVector getrow( int r ) const;
        void setrow( int r, const FloatVector& row );
        FloatVector getcolumn( int c ) const;
        void setcolumn( int c, const FloatVector& column );

        /* Generate standard matrices */
        
        void zeromatrix();
        void randommatrix();
        void unitmatrix();
        void indexmapping( const IndexMap& );
        
        /* Basic manipulation */
        
        void scale( Float );
        void set( Float );
        
        /* Special operations */
        
        DenseMatrix submatrix( const IndexMap& rows, const IndexMap& columns ) const;
        Float determinant() const;
        DenseMatrix adjunctMatrix() const;
        DenseMatrix transpose() const;
        
        /* Arithmetic operations */
        
        void add( const DenseMatrix& );
        void add( Float, const DenseMatrix& );
        void add( Float, Float, const DenseMatrix& );
        static DenseMatrix MatrixMult( const DenseMatrix&, const DenseMatrix& );
        
    private:
	
        std::vector<Float> entries;
	
};
  
  

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
    temp = DenseMatrix::MatrixMult( left, right );
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