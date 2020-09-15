#ifndef INCLUDEGUARD_DENSE_DENSEMATRIX_HPP
#define INCLUDEGUARD_DENSE_DENSEMATRIX_HPP

#include <memory>
#include <utility>
#include <vector>

class DenseMatrix;
class SparseMatrix;


#include "../basic.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"
#include "../sparse/sparsematrix.hpp"



/************************
****
****  Class for Dense Matrices 
****  - instantiates LinearOperator
****  - only basic linear arithmetics, and also matrix multiplication 
****  
************************/

class DenseMatrix final
: public LinearOperator /* every matrix is a linear operator */
{

    public:
        
        DenseMatrix( const DenseMatrix& );
        DenseMatrix( DenseMatrix&& );
        DenseMatrix& operator=( const DenseMatrix& );
        DenseMatrix& operator=( DenseMatrix&& );
        
        explicit DenseMatrix( int dim, Float initialvalue = notanumber );
        DenseMatrix( int dim, const std::function<Float(int,int)>& generator );
        DenseMatrix( int dim, const std::vector<FloatVector>& coldata );
        
        DenseMatrix( int rows, int columns, Float initialvalue = notanumber );
        DenseMatrix( int rows, int columns, const std::function<Float(int,int)>& generator );
        DenseMatrix( int rows, int columns, const std::vector<FloatVector>& coldata );
        
        explicit DenseMatrix( const ScalingOperator& );
        explicit DenseMatrix( const DiagonalOperator& );
        explicit DenseMatrix( const SparseMatrix& );
        explicit DenseMatrix( const FloatVector& );
        
        virtual ~DenseMatrix();
        
        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override 
        {
            std::shared_ptr<DenseMatrix> cloned = std::make_shared<DenseMatrix>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override 
        {
            std::unique_ptr<DenseMatrix> heir = std::make_unique<DenseMatrix>( std::move(*this) );
            return heir;
        }
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;
        
        
        DenseMatrix clone() const;
        
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;
        
        /* matrix level point of view */
        
        int numrows() const;
        int numcolumns() const;
        
        /* Access entries */
        
        Float get(int,int) const;
        void set(int,int,Float);
        Float& at( int, int ) &;
        const Float& at( int, int ) const &;
        Float& operator()( int, int ) &;
        const Float& operator()( int, int ) const &;
        
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
        
        /* Produce sparse matrix entries */
        
        // std::vector<SparseMatrix::MatrixEntry> getSparseMatrixEntries( bool clean = false ) const;
        
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
        
        DenseMatrix symmetricPart() const;
        
        DenseMatrix antisymmetricPart() const;
        
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
        
        
        bool issmall( Float eps = 0.000001 ) const;
        
        
        Float* raw();
        const Float* raw() const;

    private:
        
        Float* entries;
        
};
  




inline DenseMatrix IdentityMatrix( int dim )
{
    return DenseMatrix( dim, []( int r, int c ) -> Float{ return r==c ? 1. : 0.; } );
}




inline DenseMatrix MatrixMult( const DenseMatrix& left, const DenseMatrix& right )
{
    left.check();
    right.check();

    const int lin = left.getdimin();
    const int lout = left.getdimout();
    const int rin = right.getdimin();
    const int rout = right.getdimout();

    assert( lin == rout );

    DenseMatrix ret( lout, rin, 0. );
    
    for( int lo = 0; lo < lout; lo++ )
    for( int ri = 0; ri < rin; ri++ )
    for( int m = 0; m < rout; m++ )
        ret( lo, ri ) += left( lo, m ) * right( m, ri );

    ret.check();
    return ret;
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
    // DenseMatrix ret( left );
    // ret *= right;
    // return ret;
    return MatrixMult( left, right );
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


// inline Float weightedproduct( const DenseMatrix& W, FloatVector left, FloatVector right )
// {
//     assert( W.getdimin() == right.getdimension() );
//     assert( W.getdimout() == left.getdimension() );
//     std::vector<Float> vec( W.getdimin() * W.getdimout() );
//     for( int r = 0; r < W.getdimout(); r++ )
//     for( int c = 0; c < W.getdimin();  c++ )
//         vec[ r * W.getdimin() + c ] = left[r] * right[c] * W(r,c);
// 
//     for( int i = 0; i < vec.size(); i++ )
//     for( int j = 0; j < vec.size(); j++ )
//         if( absolute(vec[i]) > absolute(vec[j]) )
//             std::swap( vec[i], vec[j] );
// 
//     Float ret = 0.;
//     for( int i = 0; i < vec.size(); i++ )
//         ret = ret + vec[i];
// 
//     return ret;
// }




#endif
