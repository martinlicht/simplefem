#ifndef INCLUDEGUARD_SPARSE_MATRIX
#define INCLUDEGUARD_SPARSE_MATRIX

#include <vector>

class SparseMatrix;
// class MatrixEntry;

#include "../basic.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/diagonaloperator.hpp"
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

        static bool compareMatrixEntry(const MatrixEntry& first, const MatrixEntry& second) 
        {
            if( first.row < second.row )
                return true;
            else if( first.row == second.row && first.column < second.column )
                return true;
            else 
                return false;
        }

        explicit SparseMatrix( int dimout, int dimin, int numentries = 0, 
                               std::function<MatrixEntry(int)> generator = [](int i __attribute__((unused)) )->MatrixEntry{ return {0,0,0.}; } ); 
        explicit SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& );
        explicit SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& );
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& );
        virtual ~SparseMatrix();

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
        const std::vector<MatrixEntry>& getentries() const;
        void clearentries();

        int getnumberofentries() const;
        void sortentries() const;
        void sortandcompressentries() const;

        SparseMatrix getTranspose() const;
        
        virtual FloatVector apply( const FloatVector& add, Float scaling ) const override;

    private:

        mutable std::vector<MatrixEntry> entries; 
    
};



SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right );



  

#endif
