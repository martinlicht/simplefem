#ifndef INCLUDEGUARD_SPARSE_MATRIX
#define INCLUDEGUARD_SPARSE_MATRIX

#include <memory>
#include <utility>
#include <vector>

// class SparseMatrix;
// class MatrixEntry;

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"
#include "../dense/densematrix.hpp"



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

        // static bool compareMatrixEntry(const MatrixEntry& first, const MatrixEntry& second) 
        // {
        //     if( first.row < second.row )
        //         return true;
        //     else if( first.row == second.row && first.column < second.column )
        //         return true;
        //     else 
        //         return false;
        // }

        explicit SparseMatrix( int dimout, int dimin, int numentries = 0, 
                               std::function<MatrixEntry(int)> generator = [](int i __attribute__((unused)) )->MatrixEntry{ return {0,0,notanumber}; } ); 
        explicit SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& );
        explicit SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& );
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& );
        virtual ~SparseMatrix();

        SparseMatrix( const SparseMatrix& );
        SparseMatrix& operator=( const SparseMatrix& );
        SparseMatrix( SparseMatrix&& );
        SparseMatrix& operator=( SparseMatrix&& );

        virtual std::shared_ptr<LinearOperator> get_shared_pointer_to_clone() const& override {
            std::shared_ptr<SparseMatrix> cloned = std::make_shared<SparseMatrix>( *this );
            return cloned;
        }
        
        virtual std::unique_ptr<LinearOperator> get_unique_pointer_to_heir() && override {
            std::unique_ptr<SparseMatrix> heir = std::make_unique<SparseMatrix>( std::move(*this) );
            return heir;
        }
        
        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        void reserve( int ) const;
        
        const MatrixEntry& getentry( int ) const;
        MatrixEntry& getentry( int );
        void setentry( int, int, int, Float );
        void setentry( int, MatrixEntry );
        void addentry( int, int, Float );
        void addentry( MatrixEntry );
        const std::vector<MatrixEntry>& getentries() const &;
        void clearentries();

        int getnumberofentries() const;
        
        bool is_sorted( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        void sortentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        void sortandcompressentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;

        SparseMatrix getTranspose() const;
        
        FloatVector InverseDiagonalPreconditioner() const;
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;

    private:

        mutable std::vector<MatrixEntry> entries; 
    
};



SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right );

SparseMatrix SparseMatrixMultiplication( const SparseMatrix& left, const SparseMatrix& right );



inline static DiagonalOperator InverseDiagonalPreconditioner( const SparseMatrix& mat )
{

    assert( mat.getdimin() == mat.getdimout() );

    FloatVector diag( mat.getdimin(), 0. );

    const std::vector<SparseMatrix::MatrixEntry>& entries = mat.getentries();

    for( const auto& entry : entries )
        if( entry.row == entry.column ) 
            diag.at( entry.row ) += absolute( entry.value );

//     for( int i = 0; i < diag.getdimension(); i++ )
//         assert( diag.at( i ) > 0. );
    
    for( int i = 0; i < diag.getdimension(); i++ )
        if( not issmall( diag.at( i ) ) )
            diag.at( i ) = 1. / diag.at( i );
        else
            diag.at( i ) = 0.;
    
    return DiagonalOperator( mat.getdimin(), diag );

}


  

#endif
