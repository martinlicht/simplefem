#ifndef INCLUDEGUARD_SPARSE_SPARSEMATRIX_HPP
#define INCLUDEGUARD_SPARSE_SPARSEMATRIX_HPP

#include <memory>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"

class DenseMatrix;


/************************
****
****  Class for Sparse Matrices  
****  - instantiates LinearOperator
****  
************************/




class SparseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        struct MatrixEntry
        {
            int row;
            int column;
            Float value;
        };

        enum class MatrixEntrySorting {
            rowwise,
            columnwise
        };

        /* Constructors */
        
        explicit SparseMatrix( int dimout, int dimin, int numentries = 0, 
                               std::function<MatrixEntry(int)> generator = [](int i __attribute__((unused)) )->MatrixEntry{ return {0,0,notanumber}; } ); 
        // explicit SparseMatrix( int dimout, int dimin );
        explicit SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& );
        explicit SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& );
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& );
        
        /* standard interface */ 
        
        SparseMatrix() = delete;
        SparseMatrix( const SparseMatrix& );
        SparseMatrix& operator=( const SparseMatrix& );
        SparseMatrix( SparseMatrix&& );
        SparseMatrix& operator=( SparseMatrix&& );
        virtual ~SparseMatrix();

        /* standard methods for operators */
        
        virtual void check() const override;
        virtual std::string text() const override;
        virtual void printplain( std::ostream& ) const;


        /* OTHER METHODS */
        
        virtual SparseMatrix* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;

        
        /* manipulation and information */
        
        void scale ( Float s );

        bool isfinite() const;
        
        FloatVector diagonal() const;
        
        int getnumberofzeroentries() const;
        
        
        /* access and information to internal data */
        
        const std::vector<MatrixEntry>& getentries() const;
        
        std::vector<MatrixEntry>& getentries();
        
        int getnumberofentries() const;
        
        
        /* sorting entries */
        
        bool is_sorted( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        const SparseMatrix& sortentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        const SparseMatrix& sortandcompressentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;

        /* specific entry manipulations */
        
        void reserve( int ) const;
                
        const MatrixEntry& getentry( int ) const;
        
        MatrixEntry& getentry( int );
        
        void setentry( int, int, int, Float );
        
        void setentry( int, MatrixEntry );
        
        void addentry( int, int, Float );
        
        void addentry( MatrixEntry );
        
        void clearentries();
        
        /* obtain a transpose */

        SparseMatrix getTranspose() const;
        
    private:

        mutable std::vector<MatrixEntry> entries; 
    
};



SparseMatrix SparseMatrixMultiplication( const SparseMatrix& left, const SparseMatrix& right );

// FloatVector InverseDiagonalPreconditioner( const SparseMatrix& mat );

DiagonalOperator InverseDiagonalPreconditioner( const SparseMatrix& mat );



inline SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right )
{
    return SparseMatrixMultiplication( left, right );
}

inline SparseMatrix operator*( const SparseMatrix& mat, Float s )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}

inline SparseMatrix operator*( Float s, const SparseMatrix& mat )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}










  

#endif
