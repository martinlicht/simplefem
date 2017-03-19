#ifndef INCLUDEGUARD_SPARSEMATRIX
#define INCLUDEGUARD_SPARSEMATRIX

#include <list>

class SparseMatrix;
// class MatrixEntry;

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "diagonaloperator.hpp"
#include "densematrix.hpp"



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

        explicit SparseMatrix(int,int);
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& );
        virtual ~SparseMatrix();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        void addentry( int, int, Float );
        void addentry( MatrixEntry );
        const std::list<MatrixEntry>& getentries() const;
        void clearentries();

        int getnumberofentries() const;
        void sortentries() const;
        void compressentries() const;

        virtual FloatVector apply( const FloatVector& add, Float scaling ) const override;

    private:

        std::list<MatrixEntry> entries;
    
};
  
  

#endif