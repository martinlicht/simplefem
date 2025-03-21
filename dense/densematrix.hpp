#ifndef INCLUDEGUARD_DENSE_DENSEMATRIX_HPP
#define INCLUDEGUARD_DENSE_DENSEMATRIX_HPP

#include <initializer_list>
#include <utility>
#include <vector>

class DenseMatrix;
class SparseMatrix;


#include "../base/include.hpp"
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
        
        /* Constructors */
        
        explicit DenseMatrix( int dim, Float initialvalue = notanumber );
        DenseMatrix( int dim, const std::function<Float(int,int)>& generator );
        DenseMatrix( int dim, const std::vector<FloatVector>& coldata );
        DenseMatrix( int dim, const std::initializer_list<Float>& rowdata );
        
        DenseMatrix( int rows, int columns, Float initialvalue = notanumber );
        DenseMatrix( int rows, int columns, const std::function<Float(int,int)>& generator );
        DenseMatrix( int rows, int columns, const std::vector<FloatVector>& coldata );
        DenseMatrix( int rows, int columns, const std::initializer_list<Float>& rowdata );

        explicit DenseMatrix( const ScalingOperator& );
        explicit DenseMatrix( const DiagonalOperator& );
        explicit DenseMatrix( const SparseMatrix& );
        explicit DenseMatrix( const FloatVector& );
                
        DenseMatrix( int number_of_blocks, const DenseMatrix& mat, Float scaling );

        explicit DenseMatrix( const DenseMatrix&, Float scaling );
        explicit DenseMatrix( DenseMatrix&&, Float scaling );

        
        /* standard interface */ 
        
        DenseMatrix() = delete;
        DenseMatrix( const DenseMatrix& );
        DenseMatrix( DenseMatrix&& ) noexcept;
        DenseMatrix& operator=( const DenseMatrix& );
        DenseMatrix& operator=( DenseMatrix&& ) noexcept;
        virtual ~DenseMatrix() noexcept;
        
        /* standard methods for operators */

        virtual void check() const override;
        virtual std::string text() const override;
        
        std::string data_as_text( bool indexed, bool print_as_list = false ) const;
        
        /* OTHER METHODS */

        virtual DenseMatrix* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        
        
        DenseMatrix clone() const;
        
        
        using LinearOperator::apply; // import any 'apply' into the derived class' methods
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;
        
        /* matrix level point of view */
        
        int numrows() const;
        int numcolumns() const;
        
        /* Access entries */
        
        HOTCALL Float get( int r, int c ) const;
        HOTCALL void set( int r,int c, Float v );
        HOTCALL Float& at( int r, int c ) &;
        HOTCALL const Float& at( int r, int c ) const &;
        HOTCALL Float& operator()( int r, int c ) &;
        HOTCALL const Float& operator()( int r, int c ) const &;
        
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
        
        void zero_matrix();
        void random_matrix();
        void random_integer_matrix( int min, int max );
        void random_orthogonal_matrix();
        void identity_matrix();
        void indexmapping( const IndexMap& );
        
        /* Basic manipulation */
        
        void scale( Float s );
        void set( Float s );
        void add( Float s );
        
        /* Special operations */
        
        DenseMatrix submatrix( const IndexMap& rows, const IndexMap& columns ) const;
        
        /* Arithmetic operations */
        
        void add( const DenseMatrix& addendum );
        void add( Float scaling, const DenseMatrix& addendum );
        void add( Float s, Float scaling, const DenseMatrix& addendum );
        
        
        /* Simple calculations */
        
        Float min() const;
        
        Float average() const;
        
        Float sum() const;
        
        Float maxabsoluteentry() const;
        
        Float norm() const;
        
        Float maxnorm() const;
        
        Float sumnorm() const;
        
        Float l2norm() const;

        Float lpnorm( Float ) const;

        Float frobeniusnorm() const;
        
        Float norm_row_col( Float p, Float q ) const;

        Float norm_col_row( Float p, Float q ) const;

        Float norm_operator_l1() const;

        Float norm_operator_max() const;
        

        /* Measurements */
        
        FloatVector getDiagonal() const;
        
        DenseMatrix symmetricPart() const;
        
        DenseMatrix antisymmetricPart() const;
        
        Float trace() const;

        
        /* Gerschgorin circles and eigenvalue estimates */

        DenseMatrix Gerschgorin() const;
        
        DenseMatrix GerschgorinRow() const;
        
        DenseMatrix GerschgorinColumn() const;

        Float eigenvalue_estimate() const;
        
        Float operator_norm_estimate( int sample_numbers = 5, int iteration_numbers = 10 ) const;
        
        
        /* Investigations */
        
        bool is_finite() const;
        
        bool is_zero() const;
        
        bool is_positive() const;
        
        bool is_negative() const;
        
        bool is_nonnegative() const;
        
        bool is_nonpositive() const;
        
        bool is_square() const;
        
        bool is_symmetric() const;
        
        bool is_antisymmetric() const;
        
        bool is_diagonal() const;
        
        bool is_lower_left_triangular() const;
        
        bool is_lower_right_triangular() const;
        
        bool is_upper_left_triangular() const;
        
        bool is_upper_right_triangular() const;
        
        
        bool is_numerically_small( Float threshold = desired_closeness ) const;
        
        bool is_numerically_identity( Float threshold = desired_closeness ) const;
        
        
        /* Raw access */

        Float* raw();
        const Float* raw() const;

        /* Memory size */
        
        std::size_t memorysize() const;

    private:
        
        Float* entries;
        
};
  




DenseMatrix IdentityMatrix( int dim );

DenseMatrix MatrixMult( const DenseMatrix& left, const DenseMatrix& right );

DenseMatrix Conjugation( const DenseMatrix& A, const DenseMatrix& B );


DenseMatrix HilbertMatrix( int n );

DenseMatrix InvHilbertMatrix( int n );

Float HilbertDeterminant( int n );




















inline DenseMatrix operator+( const DenseMatrix& mat )
{
    return mat;
}

inline DenseMatrix operator-( const DenseMatrix& mat )
{
    DenseMatrix ret( mat, -1. ); 
    return mat;
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
    DenseMatrix temp = MatrixMult( left, right );
    left = std::move(temp);
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

inline DenseMatrix operator*( const DenseMatrix& mat, Float s )
{
    return s * mat;
}

inline DenseMatrix operator/( const DenseMatrix& mat, Float s )
{
    DenseMatrix ret( mat );
    ret /= s;
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
