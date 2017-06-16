#ifndef INCLUDEGUARD_SPARSEMATRIX
#define INCLUDEGUARD_SPARSEMATRIX

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

        explicit SparseMatrix( int dimout, int dimin );
        explicit SparseMatrix( int dimout, int dimin, int numentries );
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& );
        explicit SparseMatrix( int dimout, int dimin, int numentries, std::function<MatrixEntry(int)> generator ); 
        virtual ~SparseMatrix();

        virtual void check() const override;
        virtual void print( std::ostream& ) const override;
        virtual void printplain( std::ostream& ) const;

        const MatrixEntry& getentry( int ) const;
        MatrixEntry& getentry( int );
        void addentry( int, int, Float );
        void addentry( MatrixEntry );
        const std::vector<MatrixEntry>& getentries() const;
        void clearentries();

        int getnumberofentries() const;
        void sortentries() const;
        void compressentries() const;

        SparseMatrix getTranspose() const;
        
        virtual FloatVector apply( const FloatVector& add, Float scaling ) const override;

    private:

        std::vector<MatrixEntry> entries; // TODO: Change from list to vector 
    
};
  
  

#endif